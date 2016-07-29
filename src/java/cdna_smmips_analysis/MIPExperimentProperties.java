/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package cdna_smmips_analysis;

/**
 *
 * @author Z009157
 */
public class MIPExperimentProperties {
    Integer umiLengthExtension; // length of unique molecular identifier before extension probe

    public Integer getUmiLengthExtension() {
        return umiLengthExtension;
    }

    public Integer getUmiLengthLigation() {
        return umiLengthLigation;
    }
    Integer umiLengthLigation; // length of unique molecular identifier before extension probe
    Integer seedSequenceHashLength; // length of seed sequence to find candidate probe matches
    public MIPExperimentProperties(Integer _umiLengthLigation, Integer _umiLengthExtension, Integer _seedSequenceHashLength) {
        umiLengthExtension = _umiLengthExtension;
        umiLengthLigation = _umiLengthLigation;
        seedSequenceHashLength = _seedSequenceHashLength;
    }
}
