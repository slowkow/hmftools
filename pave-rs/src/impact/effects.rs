//! hmf-common VariantEffect and CodingEffect, and the transcript region and coding types.

/// VariantEffect, in its declaration order
#[derive(Clone, Copy, PartialEq, Eq, Debug, Hash)]
pub enum VariantEffect {
    StopGained,
    StopLost,
    StartLost,
    Frameshift,
    SpliceAcceptor,
    SpliceDonor,
    InframeInsertion,
    InframeDeletion,
    Missense,
    PhasedInframeInsertion,
    PhasedInframeDeletion,
    PhasedMissense,
    Synonymous,
    PhasedSynonymous,
    Intronic,
    FivePrimeUtr,
    ThreePrimeUtr,
    UpstreamGene,
    NonCodingTranscript,
    Other,
}

use VariantEffect::*;

pub const ALL_EFFECTS: [VariantEffect; 20] = [
    StopGained,
    StopLost,
    StartLost,
    Frameshift,
    SpliceAcceptor,
    SpliceDonor,
    InframeInsertion,
    InframeDeletion,
    Missense,
    PhasedInframeInsertion,
    PhasedInframeDeletion,
    PhasedMissense,
    Synonymous,
    PhasedSynonymous,
    Intronic,
    FivePrimeUtr,
    ThreePrimeUtr,
    UpstreamGene,
    NonCodingTranscript,
    Other,
];

/// VariantEffect.VARIANT_EFFECTS_DELIM
pub const VARIANT_EFFECTS_DELIM: &str = "&";

impl VariantEffect {
    pub fn effect(self) -> &'static str {
        match self {
            StopGained => "stop_gained",
            StopLost => "stop_lost",
            StartLost => "start_lost",
            Frameshift => "frameshift_variant",
            SpliceAcceptor => "splice_acceptor_variant",
            SpliceDonor => "splice_donor_variant",
            InframeInsertion => "inframe_insertion",
            InframeDeletion => "inframe_deletion",
            Missense => "missense_variant",
            PhasedInframeInsertion => "phased_inframe_insertion",
            PhasedInframeDeletion => "phased_inframe_deletion",
            PhasedMissense => "phased_missense",
            Synonymous => "synonymous_variant",
            PhasedSynonymous => "phased_synonymous",
            Intronic => "intron_variant",
            FivePrimeUtr => "5_prime_UTR_variant",
            ThreePrimeUtr => "3_prime_UTR_variant",
            UpstreamGene => "upstream_gene_variant",
            NonCodingTranscript => "non_coding_transcript_exon_variant",
            Other => "other",
        }
    }

    pub fn rank(self) -> i32 {
        match self {
            StopGained | StopLost | StartLost | Frameshift => 60,
            SpliceAcceptor | SpliceDonor => 50,
            InframeInsertion | InframeDeletion | Missense | PhasedInframeInsertion | PhasedInframeDeletion | PhasedMissense => 40,
            Synonymous | PhasedSynonymous => 30,
            Intronic | FivePrimeUtr | ThreePrimeUtr | UpstreamGene | NonCodingTranscript => 10,
            Other => 0,
        }
    }

    /// fromEffect: OTHER for an unknown name
    pub fn from_effect(name: &str) -> VariantEffect {
        ALL_EFFECTS.iter().copied().find(|e| e.effect() == name).unwrap_or(Other)
    }
}

pub fn is_splice(effect: VariantEffect) -> bool {
    effect == SpliceAcceptor || effect == SpliceDonor
}

pub fn is_inframe(effect: VariantEffect) -> bool {
    matches!(effect, InframeDeletion | InframeInsertion | PhasedInframeDeletion | PhasedInframeInsertion)
}

pub fn is_nonsense_or_frameshift(effect: VariantEffect) -> bool {
    matches!(effect, Frameshift | StopGained | StopLost | StartLost)
}

/// effectsToList: split at '&'
pub fn effects_to_list(effects: &str) -> Vec<VariantEffect> {
    effects.split(VARIANT_EFFECTS_DELIM).map(VariantEffect::from_effect).collect()
}

/// effectsToString
pub fn effects_to_string(effects: &[VariantEffect], delim: &str) -> String {
    effects.iter().map(|e| e.effect()).collect::<Vec<_>>().join(delim)
}

/// CodingEffect
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[allow(dead_code)]
pub enum CodingEffect {
    NonsenseOrFrameshift,
    Splice,
    Missense,
    Synonymous,
    None,
    Undefined,
}

impl CodingEffect {
    pub fn name(self) -> &'static str {
        match self {
            CodingEffect::NonsenseOrFrameshift => "NONSENSE_OR_FRAMESHIFT",
            CodingEffect::Splice => "SPLICE",
            CodingEffect::Missense => "MISSENSE",
            CodingEffect::Synonymous => "SYNONYMOUS",
            CodingEffect::None => "NONE",
            CodingEffect::Undefined => "UNDEFINED",
        }
    }

    /// CodingEffect.effect(VariantEffect)
    pub fn effect(e: VariantEffect) -> CodingEffect {
        match e {
            Frameshift | StopGained | StopLost | StartLost => CodingEffect::NonsenseOrFrameshift,
            SpliceAcceptor | SpliceDonor => CodingEffect::Splice,
            Missense | InframeDeletion | InframeInsertion | PhasedInframeDeletion | PhasedInframeInsertion | PhasedMissense => CodingEffect::Missense,
            Synonymous | PhasedSynonymous => CodingEffect::Synonymous,
            _ => CodingEffect::None,
        }
    }
}

/// TranscriptRegionType (as far as PAVE sets it)
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum RegionType {
    Unknown,
    Upstream,
    Exonic,
    Intronic,
}

impl RegionType {
    pub fn name(self) -> &'static str {
        match self {
            RegionType::Unknown => "UNKNOWN",
            RegionType::Upstream => "UPSTREAM",
            RegionType::Exonic => "EXONIC",
            RegionType::Intronic => "INTRONIC",
        }
    }
}

/// TranscriptCodingType
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum CodingType {
    Unknown,
    Coding,
    Utr5p,
    Utr3p,
    NonCoding,
    Enhancer,
}

impl CodingType {
    pub fn name(self) -> &'static str {
        match self {
            CodingType::Unknown => "UNKNOWN",
            CodingType::Coding => "CODING",
            CodingType::Utr5p => "UTR_5P",
            CodingType::Utr3p => "UTR_3P",
            CodingType::NonCoding => "NON_CODING",
            CodingType::Enhancer => "ENHANCER",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_and_ranks() {
        for e in ALL_EFFECTS {
            assert_eq!(VariantEffect::from_effect(e.effect()), e);
        }
        assert_eq!(VariantEffect::from_effect("x"), Other);
        assert_eq!(effects_to_string(&[Missense, SpliceDonor], "&"), "missense_variant&splice_donor_variant");
        assert_eq!(effects_to_list("missense_variant&splice_donor_variant"), vec![Missense, SpliceDonor]);
        assert!(Frameshift.rank() > SpliceDonor.rank());
        assert_eq!(CodingEffect::effect(PhasedInframeDeletion), CodingEffect::Missense);
    }
}
