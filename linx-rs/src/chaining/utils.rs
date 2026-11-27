//! ChainUtils and LinkFinder: assembled links, deletion bridges, possible links and chain edits.

use hmf_common::java_hash::JavaHashMap;
use hmf_common::sv::StructuralVariantType::*;

use crate::chaining::jcn_limits::jcn_match;
use crate::types::*;
use crate::utils::format_jcn;

pub const CHAIN_LINK_COUNT: usize = 0;
pub const CHAIN_ASSEMBLY_LINK_COUNT: usize = 1;
pub const CHAIN_LENGTH: usize = 2;

/// ChainUtils.breakendsAreChained: link count, assembled link count and length between two breakends
/// facing each other in the chain (counts zero if they do not)
pub fn breakends_are_chained(chain: &SvChain, var1: &Sv, v1_start: bool, var2: &Sv, v2_start: bool) -> [i32; 3] {
    let mut be1_faces_up = false;
    let mut be1_index: i32 = -1;
    let mut be2_faces_up = false;
    let mut be2_index: i32 = -1;
    let mut link_data = [0i32; 3];
    for (i, pair) in chain.linked_pairs().iter().enumerate() {
        let i = i as i32;
        if &pair.first() == var1 && pair.first_link_on_start() == v1_start {
            be1_index = i;
            be1_faces_up = true;
        } else if &pair.second() == var1 && pair.second_link_on_start() == v1_start {
            be1_index = i;
            be1_faces_up = false;
        }
        if &pair.first() == var2 && pair.first_link_on_start() == v2_start {
            be2_index = i;
            be2_faces_up = true;
        } else if &pair.second() == var2 && pair.second_link_on_start() == v2_start {
            be2_index = i;
            be2_faces_up = false;
        }
        if be1_index >= 0 || be2_index >= 0 {
            link_data[CHAIN_LINK_COUNT] += 1;
            link_data[CHAIN_LENGTH] += pair.base_length();
            if pair.is_assembled() {
                link_data[CHAIN_ASSEMBLY_LINK_COUNT] += 1;
            }
        }
        if be1_index >= 0
            && be2_index >= 0
            && ((be1_index <= be2_index && be1_faces_up && !be2_faces_up) || (be1_index >= be2_index && !be1_faces_up && be2_faces_up))
        {
            return link_data;
        }
    }
    link_data[CHAIN_LINK_COUNT] = 0;
    link_data[CHAIN_ASSEMBLY_LINK_COUNT] = 0;
    link_data
}

/// ChainUtils.foldbackChainOnLink
pub fn foldback_chain_on_link(chain: &SvChain, pair1: &Pair, pair2: &Pair) {
    let chain_breakend = if pair1.first() == pair2.first() && pair1.second() == pair2.second() {
        if pair1.first_breakend() == pair2.first_breakend() { pair1.first_breakend() } else { pair1.second_breakend() }
    } else if pair1.first() == pair2.second() && pair1.second() == pair2.first() {
        if pair1.first_breakend() == pair2.second_breakend() { pair1.first_breakend() } else { pair1.second_breakend() }
    } else {
        hmf_common::error!("chain({}) failed to add foldback pairs: {} and {}", chain.id(), pair1, pair2);
        return;
    };
    let connect_on_start = chain.open_breakend(true).as_ref() == Some(&chain_breakend);
    let link_count = chain.link_count();
    if connect_on_start {
        let existing = chain.linked_pairs();
        chain.add_link(pair1, true);
        chain.add_link(pair2, true);
        for pair in &existing {
            chain.add_link(pair, true);
        }
    } else {
        chain.add_link(pair1, false);
        chain.add_link(pair2, false);
        for index in (0..link_count).rev() {
            let pair = chain.linked_pair(index);
            chain.add_link(&pair, false);
        }
    }
}

/// ChainUtils.foldbackChainOnChain
pub fn foldback_chain_on_chain(chain: &SvChain, foldback_chain: &SvChain, pair1: &Pair, pair2: &Pair) {
    let chain_start = chain.open_breakend(true);
    let chain_end = chain.open_breakend(false);
    let fb_start = foldback_chain.open_breakend(true).unwrap();
    let fb_end = foldback_chain.open_breakend(false).unwrap();
    let other = if pair1.has_breakend(&fb_start) { pair1.other_breakend(&fb_start) } else { pair1.other_breakend(&fb_end) };
    let connect_on_start = chain_start.as_ref() == Some(&other);
    if (connect_on_start && !pair2.has_breakend(chain_start.as_ref().unwrap()))
        || (!connect_on_start && !pair2.has_breakend(chain_end.as_ref().unwrap()))
    {
        hmf_common::error!("chain({}) failed to add foldback pairs: {} and {}", chain.id(), pair1, pair2);
        return;
    }
    let existing = chain.linked_pairs();
    let fb_links = foldback_chain.linked_pairs();
    chain.add_link(pair1, connect_on_start);
    let connecting = if connect_on_start { chain_start.unwrap() } else { chain_end.unwrap() };
    let connecting_fb_start = foldback_chain.open_breakend(true).as_ref() == Some(&pair1.other_breakend(&connecting));
    if connecting_fb_start {
        for p in &fb_links {
            chain.add_link(p, connect_on_start);
        }
    } else {
        for p in fb_links.iter().rev() {
            chain.add_link(p, connect_on_start);
        }
    }
    chain.add_link(pair2, connect_on_start);
    if !connect_on_start {
        for p in existing.iter().rev() {
            chain.add_link(p, connect_on_start);
        }
    } else {
        for p in &existing {
            chain.add_link(p, connect_on_start);
        }
    }
}

/// ChainUtils.duplicateChainOnLink
pub fn duplicate_chain_on_link(chain: &SvChain, pair1: &Pair, pair2: &Pair) {
    let (end_pair, start_pair) = if chain.can_add_linked_pair_to_end(pair1) && chain.can_add_linked_pair_to_start(pair2) {
        (pair1, pair2)
    } else if chain.can_add_linked_pair_to_end(pair2) && chain.can_add_linked_pair_to_start(pair1) {
        (pair2, pair1)
    } else {
        return;
    };
    let link_count = chain.link_count();
    chain.add_link(end_pair, false);
    chain.add_link(start_pair, false);
    for index in 0..link_count {
        let pair = chain.linked_pair(index);
        chain.add_link(&pair, false);
    }
}

/// ChainUtils.reverseSectionOnBreakend
pub fn reverse_section_on_breakend(chain: &SvChain, breakend: &Be) -> bool {
    let n = chain.link_count();
    for i in 0..n {
        let pair = chain.linked_pair(i);
        if &pair.first_breakend() == breakend {
            let chain_start = chain.open_breakend(true).expect("chain start");
            let new_link = Pair::from_breakends(&chain_start, &pair.second_breakend());
            if new_link.position_distance() <= 0 {
                return false;
            }
            new_link.set_link_reason("RECIP_INV_RECONFIG", 0);
            let to_switch: Vec<Pair> = (0..i).map(|j| chain.linked_pair(j)).collect();
            {
                let mut lp = chain.linked_pairs_mut();
                for _ in 0..to_switch.len() + 1 {
                    lp.remove(0);
                }
            }
            chain.add_link(&new_link, true);
            for p in &to_switch {
                chain.add_link(p, true);
            }
            break;
        } else if &pair.second_breakend() == breakend {
            let chain_end = chain.open_breakend(false).expect("chain end");
            let new_link = Pair::from_breakends(&pair.first_breakend(), &chain_end);
            if new_link.position_distance() <= 0 {
                return false;
            }
            new_link.set_link_reason("RECIP_INV_RECONFIG", 0);
            let to_switch: Vec<Pair> = (i + 1..chain.link_count()).map(|j| chain.linked_pair(j)).collect();
            {
                let mut lp = chain.linked_pairs_mut();
                for _ in 0..to_switch.len() + 1 {
                    lp.remove(i);
                }
            }
            chain.add_link(&new_link, false);
            for p in to_switch.iter().rev() {
                chain.add_link(p, false);
            }
            break;
        }
    }
    true
}

/// ChainUtils.reconcileChains
pub fn reconcile_chains(chains: &mut Vec<Chain>, check_chain_splits: bool, mut next_chain_id: i32, use_chain_end_jcns: bool) {
    if chains.len() <= 1 {
        return;
    }
    let mut index1 = 0;
    let mut skipped_chain_closing = false;
    while index1 < chains.len().saturating_sub(1) {
        let chain1 = chains[index1].clone();
        let mut chains_merged = false;
        let mut index2 = index1 + 1;
        while index2 < chains.len() {
            let chain2 = chains[index2].clone();
            if chain1 == chain2 {
                index2 += 1;
                continue;
            }
            if check_chain_splits && identical_chain(&chain1, &chain2, false, false) {
                chain1.set_jcn_data(chain1.jcn() + chain2.jcn(), chain1.jcn_uncertainty());
                chains.remove(index2);
                chains_merged = true;
                break;
            }
            let mut jcn_matched = jcn_match(chain1.jcn(), chain1.jcn_uncertainty(), chain2.jcn(), chain2.jcn_uncertainty());
            if !jcn_matched && !check_chain_splits && !use_chain_end_jcns {
                index2 += 1;
                continue;
            }
            'be: for be1 in [SE_START, SE_END] {
                let c1_start = is_start(be1);
                let Some(breakend1) = chain1.open_breakend(c1_start) else { continue };
                for be2 in [SE_START, SE_END] {
                    let c2_start = is_start(be2);
                    let Some(breakend2) = chain2.open_breakend(c2_start) else { continue };
                    let could_join = breakend1 != breakend2 && breakend1.sv() == breakend2.sv();
                    if !could_join {
                        continue;
                    }
                    if use_chain_end_jcns && !jcn_matched {
                        if jcn_match(chain1.jcn(), chain1.jcn_uncertainty(), breakend1.jcn(), breakend1.jcn_uncertainty())
                            && jcn_match(chain2.jcn(), chain2.jcn_uncertainty(), breakend1.jcn(), breakend1.jcn_uncertainty())
                        {
                            jcn_matched = true;
                        } else if !check_chain_splits {
                            continue;
                        }
                    }
                    if !jcn_matched {
                        let new_chain = Chain::create(next_chain_id);
                        next_chain_id += 1;
                        let higher = if chain1.jcn() > chain2.jcn() { chain1.clone() } else { chain2.clone() };
                        let lower = if higher == chain1 { chain2.clone() } else { chain1.clone() };
                        new_chain.copy_from(&higher);
                        new_chain.set_jcn_data(higher.jcn() - lower.jcn(), higher.jcn_uncertainty());
                        hmf_common::debug!(
                            "splitting chain({}) jcn({}) vs chain({}) jcn({}) into new chain({}) JCN({})",
                            higher.id(),
                            format_jcn(higher.jcn()),
                            lower.id(),
                            format_jcn(lower.jcn()),
                            new_chain.id(),
                            format_jcn(new_chain.jcn())
                        );
                        higher.set_jcn_data(lower.jcn(), higher.jcn_uncertainty());
                        chains.push(new_chain);
                    }
                    hmf_common::debug!(
                        "merging chain({} links={}) {} to chain({} links={}) {}",
                        chain1.id(),
                        chain1.link_count(),
                        if c1_start { "start" } else { "end" },
                        chain2.id(),
                        chain2.link_count(),
                        if c2_start { "start" } else { "end" }
                    );
                    if c2_start {
                        for lp in chain2.linked_pairs() {
                            chain1.add_link(&lp, c1_start);
                        }
                    } else {
                        for lp in chain2.linked_pairs().iter().rev() {
                            chain1.add_link(lp, c1_start);
                        }
                    }
                    chains.remove(index2);
                    chains_merged = true;
                    break 'be;
                }
            }
            if chains_merged {
                break;
            }
            index2 += 1;
        }
        if !chains_merged {
            index1 += 1;
        } else if skipped_chain_closing {
            skipped_chain_closing = false;
            index1 = 0;
        }
    }
    let _ = &mut skipped_chain_closing;
}

/// ChainUtils.identicalChain
pub fn identical_chain(chain: &SvChain, other: &SvChain, allow_subsets: bool, allow_same_links: bool) -> bool {
    if chain.link_sum() != other.link_sum() {
        return false;
    }
    let links = chain.linked_pairs();
    let other_links = other.linked_pairs();
    if links.len() == other_links.len() {
        let exact = links.iter().zip(other_links.iter()).all(|(a, b)| a.first() == b.first() && a.second() == b.second());
        if exact {
            return true;
        }
        return allow_same_links && chain.same_links(other);
    }
    if !allow_subsets {
        return false;
    }
    let (longer, shorter) = if links.len() > other_links.len() { (&links, &other_links) } else { (&other_links, &links) };
    let mut j = 0;
    let mut matching_found = false;
    for pair in longer.iter() {
        if j >= shorter.len() {
            break;
        }
        let other_pair = &shorter[j];
        if pair.first() == other_pair.first() && pair.second() == other_pair.second() {
            matching_found = true;
            j += 1;
        } else if matching_found {
            return false;
        }
    }
    matching_found
}

// LinkFinder

/// LinkFinder.createAssemblyLinkedPairs
pub fn create_assembly_linked_pairs(cluster: &Cluster) -> Vec<Pair> {
    let mut linked_pairs: Vec<Pair> = Vec::new();
    if cluster.sv_count() < 2 {
        return linked_pairs;
    }
    let mut has_multi = false;
    let mut multi_breakends: Vec<Be> = Vec::new();
    for (_, list) in cluster.chr_breakend_entries() {
        if list.len() == 1 {
            continue;
        }
        for i in 0..list.len() - 1 {
            let lower = &list[i];
            if !lower.reverse_orient() {
                continue;
            }
            let lower_sv = lower.sv();
            if lower_sv.sv_type() == INS {
                continue;
            }
            for upper in list.iter().skip(i + 1) {
                if !upper.forward_orient() {
                    continue;
                }
                let upper_sv = upper.sv();
                if upper_sv == lower_sv {
                    continue;
                }
                if upper_sv.sv_type() == INS {
                    continue;
                }
                let (v1s, v2s) = (lower.uses_start(), upper.uses_start());
                if !have_linked_assemblies(lower_sv, upper_sv, v1s, v2s) {
                    continue;
                }
                if lower.is_assembled_link() && !contains(&multi_breakends, lower) {
                    multi_breakends.push(lower.clone());
                    has_multi = true;
                }
                if upper.is_assembled_link() && !contains(&multi_breakends, upper) {
                    multi_breakends.push(upper.clone());
                    has_multi = true;
                }
                let p = Pair::from_breakends(lower, upper);
                p.set_is_assembled();
                lower_sv.add_linked_pair(&p, v1s);
                upper_sv.add_linked_pair(&p, v2s);
                hmf_common::debug!("cluster({}) adding assembly linked pair({}) length({})", cluster.id(), p, p.base_length());
                linked_pairs.push(p);
            }
        }
    }
    if has_multi {
        has_multi = analyse_multi_connection_breakends(cluster, &multi_breakends, &mut linked_pairs);
        if has_multi {
            cluster.set_requires_replication();
        }
    }
    linked_pairs
}

fn analyse_multi_connection_breakends(cluster: &Cluster, multi: &[Be], linked_pairs: &mut Vec<Pair>) -> bool {
    let mut reassess = false;
    for breakend in multi {
        let links = breakend.sv().assembled_linked_pairs(breakend.uses_start());
        if links.len() < 2 {
            continue;
        }
        let mut spanning: Vec<Pair> = Vec::new();
        for pair in &links {
            if contains(&spanning, pair) {
                continue;
            }
            let other = pair.other_breakend(breakend);
            let other_sv = other.sv();
            if other_sv.sv_type() != DEL && other_sv.sv_type() != DUP {
                continue;
            }
            let other_links = other_sv.assembled_linked_pairs(!other.uses_start());
            if other_links.is_empty() {
                continue;
            }
            for other_pair in &links {
                if other_pair == pair {
                    continue;
                }
                let other1 = other_pair.other_breakend(breakend);
                if other_links.iter().any(|x| x.has_breakend(&other1)) {
                    hmf_common::debug!(
                        "cluster({}) breakends({} & {} type={}) both assemble to breakend({}), removing outer link({})",
                        cluster.id(),
                        breakend,
                        other,
                        other_sv.sv_type(),
                        other1,
                        other_pair
                    );
                    spanning.push(other_pair.clone());
                    reassess = true;
                }
            }
        }
        for pair in &spanning {
            breakend.sv().remove_linked_pair(pair, breakend.uses_start());
            let other = pair.other_breakend(breakend);
            other.sv().remove_linked_pair(pair, other.uses_start());
            remove_item(linked_pairs, pair);
        }
    }
    let mut has_multi = true;
    if reassess {
        has_multi = false;
        for pair in linked_pairs.iter() {
            for se in [SE_START, SE_END] {
                let b = pair.breakend_se(se);
                if linked_pairs.iter().filter(|x| *x != pair).any(|x| x.has_breakend(&b)) {
                    has_multi = true;
                    break;
                }
            }
        }
    }
    has_multi
}

/// LinkFinder.haveLinkedAssemblies: ESVEE's linked breakend ids
pub fn have_linked_assemblies(var1: &SvVarData, var2: &SvVarData, v1_start: bool, v2_start: bool) -> bool {
    let target = if v2_start { &var2.sv_data().vcf_id_start } else { &var2.sv_data().vcf_id_end };
    var1.assembly_links(v1_start).iter().any(|x| x == target)
}

/// LinkFinder.isPossibleLink
pub fn is_possible_link(chr1: &str, pos1: i32, orient1: i8, chr2: &str, pos2: i32, orient2: i8, min_distance: i32) -> bool {
    if chr1 != chr2 {
        return false;
    }
    is_possible_link_pos(pos1, orient1, pos2, orient2, min_distance)
}

pub fn is_possible_link_pos(pos1: i32, orient1: i8, pos2: i32, orient2: i8, min_distance: i32) -> bool {
    if orient1 == orient2 {
        return false;
    }
    let d = (pos1 - pos2).abs() + 1;
    if d < min_distance {
        return false;
    }
    (pos1 < pos2 && orient1 == -1) || (pos1 > pos2 && orient1 == 1)
}

/// LinkFinder.breakendsCannotLink
pub fn breakends_cannot_link(b1: &SvBreakend, b2: &SvBreakend) -> bool {
    let d = (b1.position() - b2.position()).abs() + 1;
    get_min_templated_insertion_length(b1, b2) > d
}

/// LinkFinder.getMinTemplatedInsertionLength
pub fn get_min_templated_insertion_length(b1: &SvBreakend, b2: &SvBreakend) -> i32 {
    if (b1.sv_type() == INF || b2.sv_type() == INF) && (b1.position() - b2.position()).abs() <= 1 {
        return 0;
    }
    let homology_lengths = (b1.homology().len() + b2.homology().len()) as i32;
    let min_anchor = b1.anchor_distance().min(b2.anchor_distance());
    (min_anchor - homology_lengths).max(MIN_TEMPLATED_INSERTION_LENGTH)
}

/// LinkFinder.findDeletionBridges
pub fn find_deletion_bridges(map: &JavaHashMap<String, Vec<Be>>) {
    for (_, list) in map.iter() {
        if list.len() < 2 {
            continue;
        }
        for i in 0..list.len() - 1 {
            let breakend = &list[i];
            let next = &list[i + 1];
            let (var1, var2) = (breakend.sv(), next.sv());
            if var1 == var2 {
                continue;
            }
            if breakend.orientation() == next.orientation() {
                continue;
            }
            if breakend.arm() != next.arm() {
                continue;
            }
            let base_distance = (next.position() - breakend.position() + 1) as i64;
            let min_ti = get_min_templated_insertion_length(breakend, next);
            if breakend.forward_orient() && next.reverse_orient() {
                if i + 2 < list.len() {
                    let next_next = &list[i + 2];
                    let next_base_distance = next_next.position() - next.position() + 1;
                    let next_min_ti = get_min_templated_insertion_length(next, next_next);
                    let opposite = next.orientation() != next_next.orientation();
                    if opposite && (next_base_distance as i64) < base_distance && (next_base_distance < next_min_ti || next_min_ti == 0) {
                        continue;
                    }
                }
                mark_deletion_bridge(Db::new(DbPair::new(breakend, next)));
            } else if base_distance < min_ti as i64 || min_ti == 0 {
                if have_linked_assemblies(var1, var2, breakend.uses_start(), next.uses_start()) {
                    continue;
                }
                mark_deletion_bridge(Db::new(DbPair::new(breakend, next)));
            }
        }
    }
}

fn mark_deletion_bridge(db: Db) {
    let existing1 = db.lower_sv().db_link(db.lower_link_on_start());
    let existing2 = db.upper_sv().db_link(db.upper_link_on_start());
    if existing1.as_ref().is_some_and(|e| e.length() < db.length()) || existing2.as_ref().is_some_and(|e| e.length() < db.length()) {
        return;
    }
    if let Some(e) = &existing1 {
        if e.lower_sv() == db.lower_sv() {
            e.upper_sv().set_db_link(None, e.upper_link_on_start());
        } else {
            e.lower_sv().set_db_link(None, e.lower_link_on_start());
        }
    }
    // Java compares the upper SV's existing link with the new pair's lower SV, as for the lower one
    if let Some(e) = &existing2 {
        if e.lower_sv() == db.lower_sv() {
            e.upper_sv().set_db_link(None, e.upper_link_on_start());
        } else {
            e.lower_sv().set_db_link(None, e.lower_link_on_start());
        }
    }
    db.lower_sv().set_db_link(Some(db.clone()), db.lower_link_on_start());
    db.upper_sv().set_db_link(Some(db.clone()), db.upper_link_on_start());
}
