//! MdTag: the matches, mismatches and deletions of an alignment's MD tag.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MdTagType {
    Match,
    Snv,
    Del,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MdTagElement {
    pub type_: MdTagType,
    pub length: i32,
    pub base: char,
}

impl MdTagElement {
    pub fn is_match(&self) -> bool {
        self.type_ == MdTagType::Match
    }
    pub fn is_indel(&self) -> bool {
        self.type_ == MdTagType::Del
    }
    pub fn is_snv(&self) -> bool {
        self.type_ == MdTagType::Snv
    }
}

pub struct MdTag {
    tag: String,
    elements: Vec<MdTagElement>,
}

pub const MATCH_BYTE: u8 = 0;
const DEL_TAG: char = '^';
const MATCH_TAG: char = '0';

impl MdTag {
    pub fn new(tag: &str) -> MdTag {
        let mut m = MdTag { tag: tag.to_string(), elements: Vec::new() };
        m.parse_tag();
        m
    }
    pub fn elements(&self) -> &Vec<MdTagElement> {
        &self.elements
    }
    pub fn base_length(&self) -> i32 {
        self.elements.iter().map(|x| x.length).sum()
    }
    pub fn sequence_length(&self) -> i32 {
        self.elements.iter().filter(|x| x.type_ != MdTagType::Del).map(|x| x.length).sum()
    }
    pub fn has_mismatches(&self) -> bool {
        self.elements.iter().any(|x| x.type_ != MdTagType::Match)
    }
    pub fn mismatch_count(&self) -> i32 {
        self.elements.iter().filter(|x| x.type_ != MdTagType::Match).count() as i32
    }

    pub fn extract_sub_sequence(&self, seq_index_start: i32, seq_index_end: i32, reverse: bool) -> Vec<u8> {
        let mut seq_index = 0;
        let elements: Vec<MdTagElement> = if reverse { self.elements.iter().rev().copied().collect() } else { self.elements.clone() };
        let length = seq_index_end - seq_index_start + 1;
        // a negative length is a NegativeArraySizeException in Java
        let mut match_array = vec![0u8; length as usize];
        let mut index = 0usize;
        let mut follows_del = false;
        for element in &elements {
            for _ in 0..element.length {
                if follows_del {
                    follows_del = false;
                    seq_index += 1;
                    continue;
                }
                if seq_index > seq_index_end {
                    break;
                }
                if seq_index >= seq_index_start {
                    match_array[index] = element.type_ as u8;
                    index += 1;
                    if element.type_ == MdTagType::Del {
                        follows_del = true;
                        break;
                    }
                }
                seq_index += 1;
            }
            if seq_index > seq_index_end {
                break;
            }
        }
        match_array
    }

    fn parse_tag(&mut self) {
        let mut length = 0;
        let mut prev_tag: Option<MdTagType> = None;
        for c in self.tag.chars() {
            let mut tag = if c == DEL_TAG {
                MdTagType::Del
            } else if matches!(c, 'A' | 'C' | 'G' | 'T' | 'N') {
                MdTagType::Snv
            } else {
                MdTagType::Match
            };
            if prev_tag == Some(MdTagType::Del) && tag == MdTagType::Snv {
                tag = MdTagType::Del;
            }
            let digit = || c.to_digit(10).unwrap_or_else(|| panic!("For input string: \"{c}\"")) as i32;
            if prev_tag.is_some() && Some(tag) != prev_tag {
                let prev = prev_tag.unwrap();
                if (prev == MdTagType::Match || prev == MdTagType::Del) && length > 0 {
                    let prev_char = if prev == MdTagType::Match { MATCH_TAG } else { DEL_TAG };
                    self.elements.push(MdTagElement { type_: prev, length, base: prev_char });
                }
                length = 0;
                if tag == MdTagType::Snv {
                    self.elements.push(MdTagElement { type_: tag, length: 1, base: c });
                } else if tag == MdTagType::Match {
                    length = digit();
                }
            } else if tag == MdTagType::Match {
                length = length * 10 + digit();
            } else if tag == MdTagType::Del {
                length += 1;
            }
            prev_tag = Some(tag);
        }
        if prev_tag == Some(MdTagType::Match) {
            self.elements.push(MdTagElement { type_: MdTagType::Match, length, base: MATCH_TAG });
        }
    }
}
