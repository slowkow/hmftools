//! Java package assembly.types.

pub mod assembly_link;
pub mod assembly_stats;
pub mod junction;
pub mod junction_assembly;
pub mod junction_group;
pub mod junction_sequence;
pub mod phase_group;
pub mod phase_set;
pub mod read_assembly_indices;
pub mod ref_side_soft_clip;
pub mod remote_region;
pub mod repeat_info;
pub mod support_read;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum SupportType {
    Junction,
    Indel,
    Discordant,
    JunctionMate,
    Extension,
}

impl SupportType {
    pub fn is_split_support(self) -> bool {
        matches!(self, SupportType::Junction | SupportType::Indel)
    }
    pub fn name(self) -> &'static str {
        match self {
            SupportType::Junction => "JUNCTION",
            SupportType::Indel => "INDEL",
            SupportType::Discordant => "DISCORDANT",
            SupportType::JunctionMate => "JUNCTION_MATE",
            SupportType::Extension => "EXTENSION",
        }
    }
}

/// Declared in order of strength: a lower ordinal overrides a higher one.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum AssemblyOutcome {
    LocalIndel,
    Linked,
    DupBranched,
    Secondary,
    SuppOnly,
    NoLink,
    Decoy,
    Unset,
}

impl AssemblyOutcome {
    pub fn ordinal(self) -> i32 {
        self as i32
    }
    pub fn name(self) -> &'static str {
        match self {
            AssemblyOutcome::LocalIndel => "LOCAL_INDEL",
            AssemblyOutcome::Linked => "LINKED",
            AssemblyOutcome::DupBranched => "DUP_BRANCHED",
            AssemblyOutcome::Secondary => "SECONDARY",
            AssemblyOutcome::SuppOnly => "SUPP_ONLY",
            AssemblyOutcome::NoLink => "NO_LINK",
            AssemblyOutcome::Decoy => "DECOY",
            AssemblyOutcome::Unset => "UNSET",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum RemoteReadType {
    JunctionMate,
    Discordant,
    Supplementary,
}

impl RemoteReadType {
    pub fn ordinal(self) -> usize {
        self as usize
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum LinkType {
    Split,
    Facing,
    Indel,
}

impl LinkType {
    pub fn name(self) -> &'static str {
        match self {
            LinkType::Split => "SPLIT",
            LinkType::Facing => "FACING",
            LinkType::Indel => "INDEL",
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum BaseMismatchType {
    Insert,
    Deletion,
    BaseMismatch,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum InsertionType {
    None,
    Line,
    Near,
    Other,
}

impl InsertionType {
    pub fn name(self) -> &'static str {
        match self {
            InsertionType::None => "NONE",
            InsertionType::Line => "LINE",
            InsertionType::Near => "NEAR",
            InsertionType::Other => "OTHER",
        }
    }
}

/// hmf-common StructuralVariantType
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum StructuralVariantType {
    Bnd,
    Del,
    Dup,
    Inf,
    Ins,
    Inv,
    Sgl,
}

impl StructuralVariantType {
    pub fn name(self) -> &'static str {
        match self {
            StructuralVariantType::Bnd => "BND",
            StructuralVariantType::Del => "DEL",
            StructuralVariantType::Dup => "DUP",
            StructuralVariantType::Inf => "INF",
            StructuralVariantType::Ins => "INS",
            StructuralVariantType::Inv => "INV",
            StructuralVariantType::Sgl => "SGL",
        }
    }
}

/// SvUtils.formSvType
pub fn form_sv_type(
    chr_start: &str, chr_end: &str, pos_start: i32, pos_end: i32, orient_start: crate::assembly::orientation::Orientation,
    orient_end: crate::assembly::orientation::Orientation, has_inserted_bases: bool,
) -> StructuralVariantType {
    if chr_start != chr_end {
        return StructuralVariantType::Bnd;
    }
    if orient_start != orient_end {
        let pos_diff = (pos_start - pos_end).abs();
        if pos_diff == 1 && has_inserted_bases {
            return StructuralVariantType::Ins;
        }
        if pos_diff == 0 {
            return StructuralVariantType::Dup;
        }
        let first_is_lower = pos_start < pos_end;
        if first_is_lower == orient_start.is_forward() { StructuralVariantType::Del } else { StructuralVariantType::Dup }
    } else {
        StructuralVariantType::Inv
    }
}
