//! A classification random forest as randomForest 4.7 stores it (`object$forest`), and its prediction
//! as `predict.randomForest(type = "prob")` makes it: each tree is walked from its root, going left when
//! the variable is at most the split value (classForest / predictClassTree in rf.c), and each tree votes
//! for the class of the terminal node it reaches; the probability of a class is its votes divided by the
//! votes cast.

use anyhow::{anyhow, bail, Context, Result};

use super::rfmt::parse_hex_float;

/// rf.c NODE_TERMINAL
const NODE_TERMINAL: i32 = -1;

#[derive(Clone, Debug, PartialEq)]
pub struct Node {
    pub status: i32,
    /// 1-based variable index (0 for a terminal node)
    pub var: i32,
    /// 1-based daughter nodes (0 for a terminal node)
    pub left: i32,
    pub right: i32,
    /// 1-based class of a terminal node
    pub pred: i32,
    pub split: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Forest {
    pub classes: Vec<String>,
    /// the variables in the order predict takes them (rownames of object$importance)
    pub variables: Vec<String>,
    pub trees: Vec<Vec<Node>>,
}

/// The forest exported from CHORD.rds by export_forest.R
pub const EMBEDDED_FOREST: &str = include_str!("../../resources/chord_forest.tsv");

impl Forest {
    /// The text form written by export_forest.R
    pub fn parse(text: &str) -> Result<Forest> {
        let mut classes = Vec::new();
        let mut variables = Vec::new();
        let mut ntree: Option<usize> = None;
        let mut trees: Vec<Vec<Node>> = Vec::new();
        let mut expected_nodes = 0usize;
        for (line_no, line) in text.lines().enumerate() {
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let f: Vec<&str> = line.split('\t').collect();
            let ctx = || format!("forest line {}: {line}", line_no + 1);
            match f[0] {
                "type" => {
                    if f.get(1) != Some(&"classification") {
                        bail!("{}: only classification forests are supported", ctx());
                    }
                }
                "classes" => classes = f[1..].iter().map(|s| s.to_string()).collect(),
                "variables" => variables = f[1..].iter().map(|s| s.to_string()).collect(),
                "ncat" => {
                    if f[1..].iter().any(|c| *c != "1") {
                        bail!("{}: categorical predictors are not supported", ctx());
                    }
                }
                "ntree" => ntree = Some(f.get(1).and_then(|v| v.parse().ok()).with_context(ctx)?),
                "tree" => {
                    if let Some(last) = trees.last() {
                        if last.len() != expected_nodes {
                            bail!("{}: the previous tree has {} nodes, not {expected_nodes}", ctx(), last.len());
                        }
                    }
                    expected_nodes = f.get(2).and_then(|v| v.parse().ok()).with_context(ctx)?;
                    trees.push(Vec::with_capacity(expected_nodes));
                }
                _ => {
                    if f.len() != 6 {
                        bail!("{}: expected 6 fields", ctx());
                    }
                    let int = |i: usize| f[i].parse::<i32>().with_context(ctx);
                    let node = Node {
                        status: int(0)?,
                        var: int(1)?,
                        left: int(2)?,
                        right: int(3)?,
                        pred: int(4)?,
                        split: parse_hex_float(f[5]).ok_or_else(|| anyhow!("{}: invalid split value", ctx()))?,
                    };
                    trees.last_mut().ok_or_else(|| anyhow!("{}: node before the first tree", ctx()))?.push(node);
                }
            }
        }
        if trees.last().is_some_and(|t| t.len() != expected_nodes) {
            bail!("forest: the last tree is incomplete");
        }
        let forest = Forest { classes, variables, trees };
        if ntree != Some(forest.trees.len()) {
            bail!("forest: {} trees where ntree is {ntree:?}", forest.trees.len());
        }
        forest.check()?;
        Ok(forest)
    }

    pub fn embedded() -> Result<Forest> {
        Forest::parse(EMBEDDED_FOREST)
    }

    /// Every walk ends at a terminal node with a valid class
    pub fn check(&self) -> Result<()> {
        if self.classes.is_empty() || self.variables.is_empty() || self.trees.is_empty() {
            bail!("forest without classes, variables or trees");
        }
        for (t, nodes) in self.trees.iter().enumerate() {
            for (k, n) in nodes.iter().enumerate() {
                if n.status == NODE_TERMINAL {
                    if n.pred < 1 || n.pred as usize > self.classes.len() {
                        bail!("tree {}: node {} predicts class {}", t + 1, k + 1, n.pred);
                    }
                } else {
                    let ok = |d: i32| d as usize > k + 1 && d as usize <= nodes.len();
                    if n.var < 1 || n.var as usize > self.variables.len() || !ok(n.left) || !ok(n.right) {
                        bail!("tree {}: node {} is not a valid split", t + 1, k + 1);
                    }
                }
            }
        }
        Ok(())
    }

    /// The votes for each class of the sample whose values (in `variables` order) are `x`
    pub fn votes(&self, x: &[f64]) -> Vec<u32> {
        let mut votes = vec![0u32; self.classes.len()];
        for nodes in &self.trees {
            let mut k = 0usize;
            while nodes[k].status != NODE_TERMINAL {
                let n = &nodes[k];
                let m = (n.var - 1) as usize;
                k = if x[m] <= n.split { n.left - 1 } else { n.right - 1 } as usize;
            }
            votes[(nodes[k].pred - 1) as usize] += 1;
        }
        votes
    }
}
