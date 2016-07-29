import pystan
import pandas as pd
import numpy as np
import os,sys
import statsmodels.api as sm
import argparse

################################################
####                                        #### 
#### stan_cdna_smmips_mt_np module          ####
####                                        ####
################################################


import pystan
import pandas as pd
import random
import numpy as np
import scipy.stats
from pystan import StanModel
import threading
from threading import Thread
try:
    import cPickle as pickle
except:
    sys.stderr.write("Warning: no cPickle support\n")
    import pickle

MAX_NUM_BIAS_VECTORS=2

class Stan_cdna_thread(threading.Thread):
    def __init__(self, stan_cdna_smmips):
        threading.Thread.__init__(self)
        self.stan_cdna_smmips = stan_cdna_smmips
    def run(self):
        self.stan_cdna_smmips.run_stan_model()


class Stan_cdna_smmips_mt(object):
    
    def __init__(self, mip_counts_samples_stats, general_options, column_names = None, compile_model = True):
        """
        Wrapper around Stan_cdna_smmips class so that each condition is run in a different thread.
        The bias vector is computed in this class by comparing differences between replicates from all conditions.
        These are then passed on to the individual threads.
        """
 
        self.column_names = {
                'probe_column' : 'unique_probe_id', 
                'experiment_column' : 'unique_sample_id',
                'condition_column' : 'individual_id',
                'count_column' : 'numObservedUMIs'}
        
        if not column_names is None:
            assert type(column_names) == dict
            for cn, val in column_names.iteritems():
                assert cn in self.column_names, "Column name not allowed/recognized!"
                self.column_names[cn] = val
        
        self.mip_counts_samples_stats = mip_counts_samples_stats.copy()
        self.arg_dict = general_options.copy() 

        # set normalization method
        if not "sample_normalization_method" in general_options:
            self.arg_dict['sample_normalization_method'] = 'deseq'
        print "Using sample normalization method",self.arg_dict['sample_normalization_method'] 
        # get conditions
        condition_column = self.column_names['condition_column']
        self.conditions =  list(set(mip_counts_samples_stats[condition_column]))
        num_conditions = len(self.conditions)

        # compile the model
        if compile_model:
            the_model = Stan_cdna_smmips.get_stan_code_multi_bias_model()
            self.sm = StanModel(model_code = the_model)
            self.arg_dict['StanModel'] = self.sm
 
        # make data structures for stan model
        self.stan_data_dict, self.stan_data_dict_helper = Stan_cdna_smmips.get_stan_dictionary(self.mip_counts_samples_stats, probe_column = self.column_names['probe_column'], experiment_column = self.column_names['experiment_column'], condition_column = self.column_names['condition_column'], count_column = self.column_names['count_column'], scale_size_factor = 1.0, sample_normalization_method = self.arg_dict['sample_normalization_method'])

        # we need to compute the scale factors using the data from all conditions simultaneously
        predefined_log_size_factor = {}
        for experiment in self.stan_data_dict_helper['experiment_list']:
            exp_index = self.stan_data_dict_helper['experiment_id_to_index'][experiment]
            predefined_log_size_factor[experiment] = self.stan_data_dict['log_size_factor'][exp_index]
        self.arg_dict['predefined_log_size_factor'] = predefined_log_size_factor
        self.arg_dict['predefined_lower_expression_bound'] = self.stan_data_dict['lower_expression_bound']


        # determine bias vectors
        if not (('model_probe_bias' in self.arg_dict) and (self.arg_dict['model_probe_bias'] == False)):
            # get shared bias vector
            # We do that here because we want to use the data for all the replicates to find the strongest bias vectors
            tc_diffs, df_concat, bias_matrix = Stan_cdna_smmips.get_mip_bias(self.stan_data_dict, self.stan_data_dict_helper)
            self.arg_dict['tc_diffs'] = tc_diffs
            self.arg_dict['df_concat'] = df_concat
            self.arg_dict['bias_matrix'] = bias_matrix
            print df_concat
        else:
            print "Not modelling probe bias"
            self.arg_dict['bias_matrix'] = np.zeros( (self.stan_data_dict['num_probes'],1))


        # make Stan_cdna_smmips instance for each condition
        
        self.condition_stan = {}
        self.condition_data = {}
        for condition in self.conditions:
            self.condition_data[condition] = mip_counts_samples_stats[mip_counts_samples_stats[condition_column] == condition].copy() 
            self.condition_stan[condition] = Stan_cdna_smmips(self.condition_data[condition], self.arg_dict, column_names = column_names)

        self.mcmc_finished = False

    
    def run_stan_model(self, do_mt = True):

        self.threads = {}

        if not do_mt:
            for condition in self.conditions:
                self.condition_stan[condition].run_stan_model()
        else:

            for condition in self.conditions:
                self.threads[condition] = Stan_cdna_thread(self.condition_stan[condition])
                print "Starting thread for condition", condition
                self.threads[condition].start()

            # wait until threads are completed
            for condition in self.conditions:
                self.threads[condition].join()
        
        self.mcmc_finished = True
    
    def get_conditions(self):
        return self.conditions[:]

    def get_differential_expression_samples(self, condition_1, condition_2):
        return Stan_cdna_smmips.estimate_delta(self.condition_stan[condition_1].get_expression_for_condition_samples(condition_1), self.condition_stan[condition_2].get_expression_for_condition_samples(condition_2), self.arg_dict['bias_matrix'])

    def get_differential_expression(self, condition_1, condition_2, get_std_dev = False):
        delta_samples_df = Stan_cdna_smmips.estimate_delta(self.condition_stan[condition_1].get_expression_for_condition_samples(condition_1), self.condition_stan[condition_2].get_expression_for_condition_samples(condition_2), self.arg_dict['bias_matrix'])    
        colname = "%s-%s" % (str(condition_1), str(condition_2))
        delta_df = pd.DataFrame({ colname : delta_samples_df.mean(axis=1)})
        if get_std_dev:
            delta_df_sd = pd.DataFrame({ colname : delta_samples_df.std(axis=1)})
            return delta_df, delta_df_sd
        else:
            return delta_df
 
    def get_expression_conditions(self):
        
        for cidx, condition in enumerate(self.conditions):
            df_c = self.condition_stan[condition].get_probe_expression()
            if cidx == 0:
                df_e = df_c.copy()
            else:
                assert (df_e.index != df_c.index).sum() == 0
                assert len(set(df_e.columns) & set(df_c.columns)) == 0
                df_e[df_c.columns] = df_c.copy()
        return df_e



    
    def save(self, fname_prefix):
        """
        Saves the samples from the pystan run for each condition in a dictionary in a pickle file.
        """
        assert self.mcmc_finished
        save_dict = {}
        save_dict['stan_result'] = {}
        for condition in self.conditions:
            save_dict['stan_result'][condition] = self.condition_stan[condition].stan_result
        
        save_dict['arg_dict'] = self.arg_dict
        save_dict['mip_counts_samples_stats'] = self.mip_counts_samples_stats
        save_dict['column_names'] = self.column_names

        with open(fname_prefix + ".stan_cdna_smmips_mt.pkl",'wb') as f:
            pickle.dump(save_dict, f)



    @staticmethod
    def load(fname_prefix, compile_model = False):

        with open(fname_prefix + ".stan_cdna_smmips_mt.pkl",'rb') as f:
            load_dict = pickle.load(f)

        mt = Stan_cdna_smmips_mt(load_dict['mip_counts_samples_stats'], load_dict['arg_dict'], column_names = load_dict['column_names'], compile_model = compile_model)
        for condition in load_dict['stan_result'].keys():
            mt.condition_stan[condition].stan_result = load_dict['stan_result'][condition] 

        return mt



 
class Stan_cdna_smmips(object):
    
    def __init__(self, mip_counts_samples_stats, arg_dict, column_names = None):
        self.stan_result = None
        self.mip_counts_samples_stats = mip_counts_samples_stats
        
        self.column_names = {
                'probe_column' : 'unique_probe_id', 
                'experiment_column' : 'unique_sample_id',
                'condition_column' : 'individual_id',
                'count_column' : 'numObservedUMIs'}
        
        if not column_names is None:
            assert type(column_names) == dict
            for cn, val in column_names.iteritems():
                assert cn in self.column_names, "Column name not allowed/recognized!"
                self.column_names[cn] = val
            print "Stan_cdna_smmips: ",self.column_names
        
        self.arg_dict = arg_dict.copy()
        self.setup()
        self.stan_result = None

   
    def setup(self):
        
        # default parameters necessary for setup

        # if bias_matrix was supplied through arg_dict in self.__init__,
        # then that is added to self.stan_data_dict here, where stan_run_model() will look for it.

        figure_prefix = "tmp"
        model_probe_bias = True
        show_cluster_map_probe_bias = False

        arg_dict = self.arg_dict
        
        if 'figure_prefix' in arg_dict:
            figure_prefix = arg_dict['figure_prefix']
        if 'model_probe_bias' in arg_dict:
            model_probe_bias = arg_dict['model_probe_bias']
        if 'show_cluster_map_probe_bias' in arg_dict:
            show_cluster_map_probe_bias = arg_dict['show_cluster_map_probe_bias']

        predefined_log_size_factor = {}
        predefined_lower_expression_bound = -50.0
        if 'predefined_log_size_factor' in arg_dict:
            predefined_log_size_factor = arg_dict['predefined_log_size_factor']
        if 'predefined_lower_expression_bound' in arg_dict:
            predefined_lower_expression_bound = arg_dict['predefined_lower_expression_bound']


        stan_data_dict, stan_data_dict_helper = Stan_cdna_smmips.get_stan_dictionary(self.mip_counts_samples_stats, probe_column = self.column_names['probe_column'], experiment_column = self.column_names['experiment_column'], condition_column = self.column_names['condition_column'], count_column = self.column_names['count_column'], scale_size_factor = 1.0, predefined_log_size_factor = predefined_log_size_factor, predefined_lower_expression_bound = predefined_lower_expression_bound)
 
        if model_probe_bias:
            if 'bias_matrix' in arg_dict:
                tc_diffs, df_concat, bias_matrix = arg_dict['tc_diffs'], arg_dict['df_concat'], arg_dict['bias_matrix']
                print "Using bias_matrix from arguments"
            else:
                tc_diffs, df_concat, bias_matrix = Stan_cdna_smmips.get_mip_bias(stan_data_dict, stan_data_dict_helper)
            # print df_concat.iloc[:MAX_NUM_BIAS_VECTORS,:]

            if show_cluster_map_probe_bias:
                sns.set_style("whitegrid")
                sns.set(style="ticks")
                f = sns.clustermap(tc_diffs,figsize = (10,15), row_cluster=True, col_cluster=True)
                f.savefig(figure_prefix+".tc_diff_clustermap.pdf")

        else:
            sys.stderr.write("NOT modelling probe bias\n")
            bias_matrix = np.zeros( (stan_data_dict['num_probes'], 1) )

        stan_data_dict['bias_matrix'] = bias_matrix
        stan_data_dict['num_bias_vectors'] = bias_matrix.shape[1]

        self.stan_data_dict = stan_data_dict
        self.stan_data_dict_helper = stan_data_dict_helper


    @staticmethod
    def make_dict_from_list( l, offset = 0 ):
        d = {}
        for idx, val in enumerate(l):
            assert not d.has_key(val)
            d[val] = idx + offset
        return d
    
    @staticmethod
    def get_stan_dictionary(mip_counts, probe_column = "uniqueID", experiment_column = "experiment", condition_column = "sample", count_column = "numObservedUMIs",scale_size_factor = 1.0, make_single_vector_data_layout = False, predefined_log_size_factor = {}, predefined_lower_expression_bound = -50.0, sample_normalization_method = ""):

        """ 
        Will make a stan data dictionary for all samples, experiments, etc in mip_counts.
        Each condition in the Stan model will be equivalent to a sample in the mip_counts data frame.
        """
        add_pseudo_count = 0 # this adds one to all the counts. It prevents very large negative values for probes with no reads without having to specify a strong prior.

        sd = {} # dictionary for stan
        hsd = {} # additional variables
        N = mip_counts.shape[0]   
        sd['N'] = N # number of data points

        
        #
        # probes. setup probe_id -> probe_index hash table
        #
        
        hsd['probe_list'] = sorted(list(set(mip_counts[probe_column])))
        hsd['probe_id_to_index'] = Stan_cdna_smmips.make_dict_from_list(hsd['probe_list'], offset = 0)
        sd['num_probes'] = len(hsd['probe_list'])
        
        if make_single_vector_data_layout:
            # map probes to indices for each data point
            sd['probe_index'] = [0] * N
            for i in range(N):
                sd['probe_index'][i] = hsd['probe_id_to_index'][ mip_counts[probe_column].iloc[i] ] + 1
         
        #
        # conditions
        #

        hsd['condition_list'] = sorted(list(set(mip_counts[condition_column])))
        hsd['condition_id_to_index'] = Stan_cdna_smmips.make_dict_from_list(hsd['condition_list'], offset = 0)
        sd['num_conditions'] = len(hsd['condition_list'])

        if make_single_vector_data_layout:
            # map conditions to indices for each data point
            sd['condition_index'] = [0] * N
            for i in range(N):
                sd['condition_index'][i] = hsd['condition_id_to_index'][ mip_counts[condition_column].iloc[i] ] + 1

        # replicate and experiment index
        hsd['experiment_list'] = []
        hsd['experiment_id_to_replicate_index'] = {}
        hsd['experiment_id_to_index'] = {}
        sd['experiment_list_condition_index'] = []
        sd['experiment_list_replicate_index'] = []
        
        condition_offset = 0
        for condition in hsd['condition_id_to_index'].keys():
            condition_index = hsd['condition_id_to_index'][condition]
            condition_experiments = sorted(list(set(mip_counts[mip_counts[condition_column] == condition][experiment_column])))
            hsd['experiment_id_to_replicate_index'][condition] = Stan_cdna_smmips.make_dict_from_list(condition_experiments, offset = 0)
            # note that we use replicate_experiment_id_to_index to make sure that the replicate index is in the same order
            # as the experiment index, except for the condition offset
            for experiment in condition_experiments:
                assert not hsd['experiment_id_to_index'].has_key(experiment)
                replicate_index = hsd['experiment_id_to_replicate_index'][condition][experiment]
                index = condition_offset + replicate_index
                assert not index in hsd['experiment_id_to_index'].values()
                hsd['experiment_id_to_index'][experiment] = index        
                hsd['experiment_list'].append(experiment)
                # NOTE PLUS 1 is STAN INDEX OFFSET
                sd['experiment_list_condition_index'].append(condition_index+1)
                sd['experiment_list_replicate_index'].append(replicate_index+1)
            condition_offset += len(condition_experiments)
        
        
        if make_single_vector_data_layout:
            sd['replicate_index'] = [0] * N
            for i in range(N):
                sd['replicate_index'][i] = hsd['experiment_id_to_replicate_index'][ mip_counts[condition_column].iloc[i] ][ mip_counts[experiment_column].iloc[i] ] + 1

        sd['num_experiments'] = len(hsd['experiment_list'])

        if make_single_vector_data_layout:
            # map experiments to indices for each data point
            sd['experiment_index'] = [0] * N
            for i in range(N):
                sd['experiment_index'][i] = hsd['experiment_id_to_index'][ mip_counts[experiment_column].iloc[i] ] + 1


          
        #
        # MIP unique counts
        #
        if make_single_vector_data_layout:
            sd['unique_counts'] = [0] * N
            for i in range(N):
                sd['unique_counts'][i] = mip_counts[count_column].iloc[i] + add_pseudo_count
        
        # also make table of counts probes x experiment
        table_counts = np.zeros( (sd['num_probes'], sd['num_experiments']), dtype='int' )
        for i in mip_counts.index:
            experiment_idx = hsd['experiment_id_to_index'][ mip_counts.loc[i, experiment_column] ] 
            probe_idx = hsd['probe_id_to_index'][ mip_counts.loc[i, probe_column] ]
            table_counts[probe_idx, experiment_idx] = mip_counts.loc[i, count_column] + add_pseudo_count 
        hsd['table_counts'] = table_counts
        sd['table_counts'] = table_counts.copy()
        

        if make_single_vector_data_layout:
            df = pd.DataFrame({'condition_index':sd['condition_index'],'experiment_index':sd['experiment_index'],'probe_index':sd['probe_index'] })
            hsd['df'] = df
        
        if False:
            #
            # Set size factor for each experiment
            #
            sd['log_size_factor'] = [0.0] * sd['num_experiments']
            tot_unique_counts = mip_counts.groupby(experiment_column,as_index = False)[count_column].sum()
            for i in tot_unique_counts.index:
                sd['log_size_factor'][ hsd['experiment_id_to_index'][ tot_unique_counts.loc[i, experiment_column] ] ] = np.log(tot_unique_counts.loc[i, count_column]) 
        else: 
            # alternative to compute size factors
            if predefined_log_size_factor == {}:  
                assert sample_normalization_method in ['sum','deseq','sum1m']
                print "Using sample relative normalization method: ",sample_normalization_method
                if sample_normalization_method == 'deseq':
                    lt = np.log(hsd['table_counts']+(1-add_pseudo_count)) # not clear how they deal with zeros in DESEQ.
                    sd['log_size_factor'] = scale_size_factor*np.median(lt - lt.mean(axis=1)[:,None], axis = 0)
                elif sample_normalization_method == 'sum':
                    log_sum = np.log(hsd['table_counts'].sum(axis=0))
                    sd['log_size_factor'] = log_sum - log_sum.mean()
                elif sample_normalization_method == 'sum1m':
                    log_sum = np.log((hsd['table_counts']-add_pseudo_count).sum(axis=0))
                    sd['log_size_factor'] = +log_sum - np.log(1e6)
                print "WARNING manual setting of lower expression bound"
                sd['lower_expression_bound'] = -10.0 # FIXME np.min([-50.0, sd['log_size_factor'].min()*10.0])
            else:
                print "Using pre-defined size factor"
                lsf = []
                for experiment in hsd['experiment_list']:
                    lsf.append( predefined_log_size_factor[experiment] )
                sd['log_size_factor'] = np.array(lsf, dtype = 'float')
                sd['lower_expression_bound'] = predefined_lower_expression_bound


        hsd['add_pseudo_count'] = add_pseudo_count 
        assert add_pseudo_count == 0
        print "Lower expression bound: %1.2g lowest log_size_factor: %1.2g"  % (sd['lower_expression_bound'], sd['log_size_factor'].min())
        
        return sd, hsd
    
    @staticmethod
    def get_mip_bias(stan_data_dict, stan_data_dict_helper, pval_threshold = 1e-5):
        # Estimates the mip probe bias vector from technical or biological replicates.
        # Correlate the difference between one pair of replicates with the 
        # difference between another pair of replicates.
        # The pair of pairs with the highest correlation is then used to estimate the
        # probe bias. Nosum(axis=0).median()
        # However, the difference is calculated only between pairs from the same condition.
        
        # For this version we require that no pseudo counts were added
        assert stan_data_dict_helper['add_pseudo_count'] == 0  
        add_pseudo_count = stan_data_dict_helper['add_pseudo_count']
        if True:
            #tc =  stan_data_dict['table_counts'].copy()
            #tc_med = np.median(tc.sum(axis=0))
            #tc = np.log2(tc/(1.0*tc.sum(axis=0))*tc_med)
            tc =  np.log2(stan_data_dict['table_counts'] + (1-add_pseudo_count) )
            tc_med = np.median(tc.sum(axis=0))
            tc = (tc/(1.0*tc.sum(axis=0))*tc_med)
        else:
            tc = np.log2(stan_data_dict['table_counts'] + (1-add_pseudo_count) )
        
        num_probes = tc.shape[0]
        
        hsd = stan_data_dict_helper
        sd = stan_data_dict
        
        # make list of lists contain index of experiment in table_counts
        replicate_list = []
        for i in range(sd['num_conditions']):
            replicate_list.append([])
            
        for exp_idx, condition_index in enumerate(sd['experiment_list_condition_index']):
            replicate_list[ condition_index - 1 ].append(exp_idx)
        print exp_idx, condition_index, replicate_list

        num_pairs = 0
        for rl in replicate_list:
            num_pairs += len(rl)*(len(rl)-1)/2

        tc_diffs = np.zeros( (num_probes, num_pairs))
        k = 0
        tc_diffs_pairs = []
        for il in replicate_list:
            for i in il:
                for j in il:
                    if j > i:
                        tc_diffs[:,k] = (tc[:,i]-tc[:,j])
                        tc_diffs_pairs.append([k,i,j])
                        k += 1
        k = 0
        df = pd.DataFrame(columns = ["p1.1","p1.2","p2.1","p2.2","corr"])
        for i in range(num_pairs):
            for j in range(i+1, num_pairs):
                if len(set(tc_diffs_pairs[i][1:])&set(tc_diffs_pairs[j][1:]))==0:
                    # only do non-overlapping pairs of pairs
                    cor = scipy.stats.spearmanr(tc_diffs[:,i], tc_diffs[:,j])
                    df.loc[k, "p1.1"] = tc_diffs_pairs[i][1]
                    df.loc[k, "p1.2"] = tc_diffs_pairs[i][2]
                    df.loc[k, "p2.1"] = tc_diffs_pairs[j][1]
                    df.loc[k, "p2.2"] = tc_diffs_pairs[j][2]
                    df.loc[k, "k1"] = tc_diffs_pairs[i][0]
                    df.loc[k, "k2"] = tc_diffs_pairs[j][0]
                    df.loc[k, "corr"] = cor[0]            
                    df.loc[k, "pval_corr"] = cor[1]            
                    k += 1
        # print df
        # greedy select pairs of pairs with high correlations 
        df = df.sort_values("pval_corr")
        df_concat = pd.DataFrame(columns=df.columns)
        df_concat.loc[0] = df.iloc[0,:]
        pairs = list(df[["p1.1","p1.2","p2.1","p2.2"]].iloc[0,:])
        i = 1
        while i<25:
            df_new = df[ ~(df["p1.1"].isin(pairs) | df["p1.2"].isin(pairs) | df["p2.1"].isin(pairs) | df["p2.2"].isin(pairs))]
            if df_new.shape[0] == 0:
                break
            df_concat.loc[i] = df_new.iloc[0,:]
            pairs = list(set(pairs) | set(df_new[["p1.1","p1.2","p2.1","p2.2"]].iloc[0,:]))
            df = df_new
            i += 1
        
        df_concat = df_concat[ df_concat["pval_corr"]<pval_threshold ]
        for c in ["k1","k2"]:
            df_concat[c] = pd.Series(df_concat[c],dtype='int')
        
        # df_concat now contains rows of correlations between pairs of differences.
        # We will make a list of bias vectors for each row, and let the MCMC model sort out the contribution from each.
        # TODO: replace with more powerful mixture modelling/clustering
        if df_concat.shape[0] == 0:
            num_bias_vectors = 1
            bias_matrix = np.zeros( (num_probes, num_bias_vectors))
        else:
            num_bias_vectors = min([df_concat.shape[0], MAX_NUM_BIAS_VECTORS]) # number of rows
            bias_matrix = np.zeros( (num_probes, num_bias_vectors))
            vidx = 0
            for i in df_concat.index[:num_bias_vectors]:
                corr = df_concat.loc[i, "corr"]
                b1 = tc_diffs[:, df_concat["k1"][i]]
                b2 = tc_diffs[:, df_concat["k2"][i]]
                if corr < 0:
                    bias_vec = (b1-b2)/2.0
                else:
                    bias_vec = (b1+b2)/2.0
                bias_matrix[:,vidx] = bias_vec
                vidx += 1
             
        return tc_diffs, df_concat, bias_matrix

    @staticmethod
    def estimate_delta(probe_expression_samples_condition_1,probe_expression_samples_condition_2, probe_bias):
        """
        Given set of MCMC samples estimates differences between condition after regressing out probe bias
        """
        assert not np.any(probe_expression_samples_condition_1.shape != probe_expression_samples_condition_2.shape)
        assert not np.any(probe_expression_samples_condition_1.index != probe_expression_samples_condition_2.index)
        S = probe_expression_samples_condition_1.shape[1]
        P = probe_expression_samples_condition_1.shape[0]
        delta = pd.DataFrame( np.zeros( (P,S) ), index = probe_expression_samples_condition_1.index )
        columns = []
        for i in range(probe_bias.shape[1]):
            columns.append("probe_bias_%d" % i)
        
        for s in range(S):
            df = pd.DataFrame({'diff':probe_expression_samples_condition_1.iloc[:,s]-probe_expression_samples_condition_2.iloc[:,s]})
            for i in range(probe_bias.shape[1]):
                df[columns[i]] = probe_bias[:,i]
            #X = sm.add_constant(df[columns]) # Only do this if you expect mean change to be zero.
            X = df[columns]
            m = sm.OLS(df['diff'], X).fit()
            delta[s] = df['diff'] - m.predict(X)
        return delta

    @staticmethod
    def make_stan_result(stan_fit, stan_data_dict, stan_data_dict_helper, arg_dict, mip_counts_samples_stats):
        df_dict = {}
        df_dict['experiment_overdispersion'] = pd.DataFrame(stan_fit['experiment_phi'], columns = stan_data_dict_helper['experiment_list'])
        df_dict['experiment_bias_coefficient'] = {'experiment_list' : stan_data_dict_helper['experiment_list'], 'samples_by_bias_vectors_by_experiments' : stan_fit['coef_bias']}
        
        # substract the threshold again
        probe_expression = stan_fit['probe_expression']+stan_data_dict['lower_expression_bound']
        
        # probes with zero counts in all experiments will have arbitrarily large negative values, which causes issues when computing log-fold change.
        # We set the expression of these probes to zero. Essentially this puts a cap on differential expression for probes who have zero expression in one
        # condition and very high expression in the other.
        probe_expression[probe_expression<0.0] = 0.0
        df_dict['condition_probe_expression'] = { 'condition_list' : stan_data_dict_helper['condition_list'], 'probe_list' : stan_data_dict_helper['probe_list'], 'samples_by_conditions_by_probes' : probe_expression}
        df_dict['stan_data_dict'] = stan_data_dict
        df_dict['stan_data_dict_helper'] = stan_data_dict_helper
        df_dict['arg_dict'] = arg_dict
        df_dict['mip_counts_samples_stats'] = mip_counts_samples_stats
        return df_dict
        
    @staticmethod
    def get_stan_code_multi_bias_model():
        the_model = """
            data { 
                int<lower=1> num_probes;      // number of smMIP probes
                int<lower=1> num_conditions;  // condition can be different individual for same cell type, or different cell line etc.
                int<lower=1> N;               // number of data points in data frame 
                int<lower=1> num_experiments; // can have multiple experiments per condition (e.g. replicates, or capture at different input concentrations) 
                int<lower=1> num_bias_vectors;
                int<lower=0> experiment_list_replicate_index[num_experiments];
                int<lower=0> experiment_list_condition_index[num_experiments];

                real         log_size_factor[num_experiments]; // size factor (e.g. log-number-of-reads for each experiment). Scaling factor for mean.
                int<lower=0> table_counts[num_probes, num_experiments]; // actual observations 
                matrix[num_probes,num_bias_vectors] bias_matrix;
                real lower_expression_bound; // lower bound on expression
             } 

             parameters { 
                 real experiment_phi[num_experiments]; // negative binomial overdispersion factor
                 vector[num_bias_vectors] coef_bias[num_experiments];
                 // matrix<lower=0>[num_conditions, num_probes] probe_expression;
                 matrix<lower=0>[num_conditions, num_probes] probe_expression;
             } 

             model { 
                 int rep_index;
                 int condition_index;
                 for (e in 1:num_experiments) {
                     experiment_phi[e] ~ gamma(1,0.10);
                     coef_bias[e] ~ normal(0.0, 10.0);
                 }

                 for (e in 1:num_experiments) {
                     rep_index <- experiment_list_replicate_index[e];
                     condition_index <- experiment_list_condition_index[e]; 
                     if (rep_index != 1) {
                         for (p in 1:num_probes) {
                             table_counts[p,e] ~ neg_binomial_2(exp(log_size_factor[ e ] + lower_expression_bound + bias_matrix[ p ]* coef_bias[ e ] + probe_expression[ condition_index, p ]), experiment_phi[e]); 
                         }
                     } else {
                         for (p in 1:num_probes) {
                             table_counts[p,e] ~ neg_binomial_2(exp(log_size_factor[ e ] + lower_expression_bound + probe_expression[ condition_index, p ]), experiment_phi[e]);                  
                        }
                     }
                 }
             } 
            """
        return the_model

    #def run_stan_model(self, seed = 222672, num_stan_iter = 1000, figure_prefix = "tmp", model_probe_bias = True, show_cluster_map_probe_bias = False):
    def run_stan_model(self):
        
        # default parameters
        seed = 222672
        num_stan_iter = 1000

        if 'seed' in self.arg_dict: 
            seed = self.arg_dict['seed'] 
        if 'num_stan_iter' in self.arg_dict:
            num_stan_iter = self.arg_dict['num_stan_iter']
        
        def myinit():
            params={}
            
            elci = self.stan_data_dict['experiment_list_condition_index']
            table_counts = 0.01 + np.log(self.stan_data_dict_helper['table_counts']+(1-self.stan_data_dict_helper['add_pseudo_count'])) - (self.stan_data_dict['lower_expression_bound'])

            conditions = set(elci)
             
            probe_expr_init = np.zeros( (self.stan_data_dict['num_probes'], self.stan_data_dict['num_conditions']) )

            condition_to_experiments = {}
            for experiment_index in range(self.stan_data_dict['num_experiments']):
                condition_index =  elci[experiment_index] # note condition_index starts at 1
                condition_to_experiments[ condition_index ] = condition_to_experiments.get( condition_index, [] )
                condition_to_experiments[ condition_index ].append(experiment_index)
            
            for condition_index in conditions:
                assert condition_index >= 1
                probe_expr_init[:, condition_index-1] = table_counts[:, condition_to_experiments[ condition_index ] ].mean(axis=1)
            
            params['probe_expression']= probe_expr_init.T
            params['experiment_phi']=np.ones( (self.stan_data_dict['num_experiments']) )*0.5
            params['coef_bias']=np.ones( (self.stan_data_dict['num_experiments'], self.stan_data_dict['num_bias_vectors']) )*0.5
            return params
        

        init = myinit()
        if not 'StanModel' in self.arg_dict:
            the_model = Stan_cdna_smmips.get_stan_code_multi_bias_model()
            self.sm = StanModel(model_code = the_model)
        else:
            self.sm = self.arg_dict['StanModel']

        
        if True:
            # tmp = self.stan_data_dict['table_counts'] 
            # self.stan_data_dict['table_counts'] = tmp - 1
            fit = self.sm.sampling(data=self.stan_data_dict, iter=num_stan_iter, chains=2, verbose = True, seed = seed, init = myinit)
            print fit
            fit_e = fit.extract()
            self.stan_result = Stan_cdna_smmips.make_stan_result(fit_e, self.stan_data_dict, self.stan_data_dict_helper, self.arg_dict, self.mip_counts_samples_stats)
        else:
            self.stan_result = None

        return self.stan_result

    def get_differential_expression(self, condition_1, condition_2):
        delta_samples_df = Stan_cdna_smmips.estimate_delta(self.get_expression_for_condition_samples(condition_1), self.get_expression_for_condition_samples(condition_2), self.stan_result['stan_data_dict']['bias_matrix'])    
        colname = "%s-%s" % (str(condition_1), str(condition_2))
        delta_df = pd.DataFrame({ colname : delta_samples_df.mean(axis=1)})
        return delta_df
 
    def get_expression_for_condition_samples(self, condition):
        """
        Returns pandas DF of probes by samples for scaled log expression for the condition
        """
        assert not self.stan_result is None, "Run the model first"
        condition_index = self.stan_result['condition_probe_expression']['condition_list'].index(condition)
        condition_expression = pd.DataFrame(self.stan_result['condition_probe_expression']['samples_by_conditions_by_probes'][:,condition_index,:].T, index = self.stan_result['condition_probe_expression']['probe_list'])
        return condition_expression

    def get_probe_expression(self):
        return pd.DataFrame(self.stan_result['condition_probe_expression']['samples_by_conditions_by_probes'].mean(axis=0).T, index = self.stan_result['condition_probe_expression']['probe_list'], columns = self.stan_result['condition_probe_expression']['condition_list'] )
    
    @staticmethod
    def save_to_pickles(pickle_fname_prefix, inst):
        model_pickle_fname = pickle_fname_prefix + ".model.pkl"
        instance_pickle_fname = pickle_fname_prefix + ".instance.pkl"

        with open(model_pickle_fname,'wb') as f:
            pickle.dump(inst.sm, f)

        with open(instance_pickle_fname,'wb') as f:
            pickle.dump(inst, f)
    
    @staticmethod
    def load_from_pickles(pickle_fname_prefix):
        model_pickle_fname = pickle_fname_prefix + ".model.pkl"
        instance_pickle_fname = pickle_fname_prefix + ".instance.pkl"
        # according to PyStan, model needs to be unpickled before unpickling the fit object
        with open(model_pickle_fname,'rb') as f:
            sm = pickle.load(f)

        with open(instance_pickle_fname,'rb') as f:
            stan_cdna_smmips = pickle.load(f)

        return stan_cdna_smmips

    def save_result_to_pickle(self, pickle_fname_prefix):
        result_pickle_fname = pickle_fname_prefix + ".result.pkl"

        with open(result_pickle_fname,'wb') as f:
            pickle.dump(self.stan_result, f)

    def load_result_from_pickle(self, pickle_fname_prefix):
        result_pickle_fname = pickle_fname_prefix + ".result.pkl"

        with open(result_pickle_fname,'rb') as f:
            self.stan_result = pickle.load(f)

        self.stan_data_dict_helper = self.stan_result['stan_data_dict_helper']
        self.stan_data_dict = self.stan_result['stan_data_dict']
        self.arg_dict = self.stan_result['arg_dict']
        self.mip_counts_samples_stats = self.stan_result['mip_counts_samples_stats']




################################################
####                                        #### 
#### stan_cdna_smmips_mt_np module          ####
####                                        ####
################################################




# helper functions
def get_sample_from_merged_df(df, sample = ""):
    dfs = df[df['unique_sample_id']==sample].copy()
    return dfs

def normalize_by_sample(df, col = "numObservedUMIsCorrected", normF = 10000.0):
    dfn = df.copy()
    ncol = col + "Norm"
    samples = list(df['unique_sample_id'].unique())
    for sample in samples:
        
        i = dfn[dfn["unique_sample_id"] == sample].index
        nu = dfn.loc[i,col].sum()
        dfn.loc[i, ncol] = df.loc[i, col]/float(nu)*float(normF)
    return dfn



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output_simple_normalization", help="Output simple log2(molecule per million) expression values", default = False)
    parser.add_argument("--compare_conditions", help ="comma-separated list of pairs of conditions for which differential expression should be estimated. For example: condition_A/condition_B,condition_A/condition_C")
    parser.add_argument("--save_model_as_pickle", help ="save Python model based on the data as pickle", action = "store_true")
    parser.add_argument("--load_model_from_pickle", help ="Do not perform MCMC inference but instead load MCMC results from previous run/pickle. (At your own risk).", action="store_true")
    parser.add_argument("--excel", help = "save output files also as excel file", action = "store_true")
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("--merged_counts_file", help="File with merged counts from all experie", required = True)
    requiredNamed.add_argument("--output_file_prefix", help ="Prefix of output files. ", required = True)

    args = parser.parse_args()

    # define output file names
    pickle_output_file = args.output_file_prefix + ".stan_model.pkl"
    bayesian_expression_file = args.output_file_prefix + ".bayesian_estimate_log_expression.txt"
    bayesian_differential_expression_file = args.output_file_prefix + ".bayesian_estimate_log2_differential_expression.txt"
    

    # counts
    print "Reading merged counts file",args.merged_counts_file
    counts = pd.read_table(args.merged_counts_file, sep = "\t")
    assert "condition" in counts.columns
    assert "experiment" in counts.columns

    conditions = list(set(counts["condition"]))
    print "Detected the following conditions in merged counts file:"," ".join(conditions)
    
    # check if all conditions for which differential expression should be tested are present in the count file
    if args.compare_conditions is not None:
        diff_condition_list = []
        diff_condition_pairs = []
        for pair in args.compare_conditions.split(','):
            pair_split = pair.split('/')
            if len(pair_split) == 1 and pair_split[0] == '':
                sys.stderr.write("Warning, empty pairs of condition\n")
            elif len(pair_split) != 2:
                sys.stderr.write("Error: the two conditions that should be compared must be separated by a '/'\n")
                sys.stderr.write("       condition pair: %s\n" % pair)
                sys.exit(1)
            else:
                diff_condition_list += pair_split
                diff_condition_pairs.append(pair_split)

        diff_condition_list = list(set(diff_condition_list))
        all_conditions_present = True
        for condition in diff_condition_list:
            if condition not in conditions:
                sys.stderr.write("Error: the condition %s specified for the --compare_conditions options was not found in the merged counts file %s\n" % (condition, args.merged_counts_file))
                all_conditions_present = False
        if not all_conditions_present:
            sys.stderr.write("Aborting. Not all conditions specified for --compare_conditions are present in merged counts file %s\n" % args.merged_counts_file)
            sys.exit(2)

    # run Bayesian model
    
    if args.load_model_from_pickle:
        print "Loading MCMC results from pickle",pickle_output_file
        stan = Stan_cdna_smmips_mt.load(pickle_output_file)
    else:
        print "Constructing Bayesian model based on the conditions in the merged counts file"
        stan = Stan_cdna_smmips_mt(counts, {'num_stan_iter':1000,
                                            'model_probe_bias': True,
                                            'sample_normalization_method': 'deseq'},
                               column_names = {'condition_column':'condition',
                                               'count_column':'numObservedUMIsCorrected',
                                               'experiment_column':'experiment'})
        print "Performing MCMC inference with the Bayesian model"
        stan.run_stan_model()
         
        if args.save_model_as_pickle:
            print "Saving Stan results and Bayesian model to", pickle_output_file
            stan.save(pickle_output_file)

    print "Output Bayesian model estimate log-expression for all conditions to", bayesian_expression_file
    bayes_log_expr = stan.get_expression_conditions()
    bayes_log_expr.to_csv(bayesian_expression_file,sep = "\t")
    if args.excel:
        bayes_log_expr.to_excel(bayesian_expression_file.replace('.txt','.xls'))

    # estimate differential expression
    if args.compare_conditions is not None:
        print "Estimating differential expression"
        for pair_idx, pair in enumerate(diff_condition_pairs):
            condition1 = pair[0]
            condition2 = pair[1]
            print "\tEstimating differential expression for condition %s vs %s" % (condition1, condition2)
            de_samples = stan.get_differential_expression_samples(condition1, condition2)
            
            # de_samples has shape (probes, MCMC samples)
            de_mean = de_samples.mean(axis=1)/np.log(2.0) # series
            de_sd = de_samples.std(axis=1)/np.log(2.0) # series
            mean_label = 'mean_log2_DE_%s_%s' % (condition1, condition2)
            sd_label = 'stdev_log2_DE_%s_%s' % (condition1, condition2)
            if pair_idx == 0:
                de_df = pd.DataFrame({mean_label : de_mean, sd_label : de_sd})
            else:
                de_df[mean_label] = de_mean
                de_df[sd_label] = de_df
            
        print "\twriting results to file",bayesian_differential_expression_file
        de_df.to_csv(bayesian_differential_expression_file, sep = "\t")
        if args.excel:
            de_df.to_excel(bayesian_differential_expression_file.replace('.txt','.xls'))
    else:
        print "No conditions to estimate differential expression for were specified."
    
    print "Finished"


if __name__ == "__main__":
    main()
