use anyhow::{ensure, Result};
use std::{cmp::Reverse, collections::BinaryHeap};

pub trait GraphIndex: Copy + Ord + std::fmt::Debug + Send + Sync {
    fn from_usize(value: usize) -> Self;
    fn as_usize(self) -> usize;
}
impl GraphIndex for usize {
    #[inline]
    fn from_usize(value: usize) -> Self {
        value
    }
    #[inline]
    fn as_usize(self) -> usize {
        self
    }
}
impl GraphIndex for u32 {
    #[inline]
    fn from_usize(value: usize) -> Self {
        debug_assert!(value <= u32::MAX as usize);
        value as u32
    }
    #[inline]
    fn as_usize(self) -> usize {
        self as usize
    }
}

#[derive(Debug)]
pub struct TypedPeelGraph<I: GraphIndex> {
    pub positions: Vec<I>,
    pub order: Vec<(I, I)>,
    pub dropped: Vec<I>,
}

/// Checked narrow graph storage, with an unrestricted native-index fallback.
pub enum PeelGraph {
    Compact(TypedPeelGraph<u32>),
    Wide(TypedPeelGraph<usize>),
}
impl PeelGraph {
    pub fn dropped_len(&self) -> usize {
        match self {
            Self::Compact(g) => g.dropped.len(),
            Self::Wide(g) => g.dropped.len(),
        }
    }
    pub fn dropped(&self) -> impl Iterator<Item = usize> + '_ {
        (0..self.dropped_len()).map(|i| match self {
            Self::Compact(g) => g.dropped[i] as usize,
            Self::Wide(g) => g.dropped[i],
        })
    }
    pub fn assign_verified(
        &mut self,
        keys: &[u64],
        values: &[u64],
        slots: usize,
        hashes: usize,
    ) -> Result<Vec<u64>> {
        match self {
            Self::Compact(g) => assign_verified(g, keys, values, slots, hashes),
            Self::Wide(g) => assign_verified(g, keys, values, slots, hashes),
        }
    }
}
fn assign_verified<I: GraphIndex>(
    graph: &mut TypedPeelGraph<I>,
    keys: &[u64],
    values: &[u64],
    slots: usize,
    hashes: usize,
) -> Result<Vec<u64>> {
    ensure!(
        keys.len() == values.len()
            && values.len().checked_mul(hashes) == Some(graph.positions.len()),
        "invalid graph dimensions"
    );
    graph.dropped.sort_unstable();
    ensure!(
        graph.dropped.windows(2).all(|p| p[0] < p[1]),
        "duplicate dropped equations"
    );
    ensure!(
        graph
            .dropped
            .last()
            .is_none_or(|i| i.as_usize() < values.len()),
        "invalid dropped equation"
    );
    let mut cells = vec![0; slots];
    for &(edge, pivot) in graph.order.iter().rev() {
        let edge = edge.as_usize();
        let pivot = pivot.as_usize();
        let mut value = values[edge];
        for &slot in &graph.positions[edge * hashes..(edge + 1) * hashes] {
            let slot = slot.as_usize();
            if slot != pivot {
                value ^= cells[slot];
            }
        }
        cells[pivot] = value;
    }
    let mut dropped = graph.dropped.iter().peekable();
    for (edge, &value) in values.iter().enumerate() {
        if dropped.peek().is_some_and(|i| i.as_usize() == edge) {
            dropped.next();
            continue;
        }
        let actual = graph.positions[edge * hashes..(edge + 1) * hashes]
            .iter()
            .fold(0, |value, slot| value ^ cells[slot.as_usize()]);
        ensure!(
            actual == value,
            "incorrect merged payload for key {}",
            keys[edge]
        );
    }
    Ok(cells)
}

/// Same observable minimum-degree / highest-cell-ID ordering as the old heap.
pub trait DegreeQueue {
    fn new(slots: usize) -> Self;
    fn push(&mut self, item: (Reverse<u32>, usize));
    fn pop(&mut self) -> Option<(Reverse<u32>, usize)>;
    fn peek(&self) -> Option<(Reverse<u32>, usize)>;
    fn update(&mut self, cell: usize, degree: u32);
}
impl DegreeQueue for BinaryHeap<(Reverse<u32>, usize)> {
    fn new(slots: usize) -> Self {
        Self::with_capacity(slots)
    }
    fn push(&mut self, item: (Reverse<u32>, usize)) {
        BinaryHeap::push(self, item);
    }
    fn pop(&mut self) -> Option<(Reverse<u32>, usize)> {
        BinaryHeap::pop(self)
    }
    fn peek(&self) -> Option<(Reverse<u32>, usize)> {
        BinaryHeap::peek(self).copied()
    }
    fn update(&mut self, cell: usize, degree: u32) {
        if degree > 1 {
            BinaryHeap::push(self, (Reverse(degree), cell));
        }
    }
}

pub struct IndexedDegreeHeap {
    nodes: Vec<(u32, u32)>, // degree, cell
    location: Vec<u32>,
}
impl IndexedDegreeHeap {
    #[inline]
    fn higher(a: (u32, u32), b: (u32, u32)) -> bool {
        a.0 < b.0 || (a.0 == b.0 && a.1 > b.1)
    }
    fn swap(&mut self, a: usize, b: usize) {
        self.nodes.swap(a, b);
        self.location[self.nodes[a].1 as usize] = a as u32;
        self.location[self.nodes[b].1 as usize] = b as u32;
    }
    fn up(&mut self, mut position: usize) {
        while position > 0 {
            let parent = (position - 1) / 2;
            if !Self::higher(self.nodes[position], self.nodes[parent]) {
                break;
            }
            self.swap(position, parent);
            position = parent;
        }
    }
    fn down(&mut self, mut position: usize) {
        loop {
            let left = position * 2 + 1;
            if left >= self.nodes.len() {
                break;
            }
            let right = left + 1;
            let child =
                if right < self.nodes.len() && Self::higher(self.nodes[right], self.nodes[left]) {
                    right
                } else {
                    left
                };
            if !Self::higher(self.nodes[child], self.nodes[position]) {
                break;
            }
            self.swap(child, position);
            position = child;
        }
    }
    fn remove(&mut self, position: usize) -> (u32, u32) {
        let item = self.nodes[position];
        let last = self.nodes.pop().unwrap();
        self.location[item.1 as usize] = u32::MAX;
        if position < self.nodes.len() {
            self.nodes[position] = last;
            self.location[last.1 as usize] = position as u32;
            if position > 0 && Self::higher(last, self.nodes[(position - 1) / 2]) {
                self.up(position);
            } else {
                self.down(position);
            }
        }
        item
    }
}
impl DegreeQueue for IndexedDegreeHeap {
    fn new(slots: usize) -> Self {
        Self {
            nodes: Vec::with_capacity(slots),
            location: vec![u32::MAX; slots],
        }
    }
    fn push(&mut self, (Reverse(degree), cell): (Reverse<u32>, usize)) {
        let position = self.location[cell];
        if position == u32::MAX {
            let position = self.nodes.len();
            self.nodes.push((degree, cell as u32));
            self.location[cell] = position as u32;
            self.up(position);
        } else {
            let position = position as usize;
            let old = self.nodes[position];
            self.nodes[position] = (degree, cell as u32);
            if Self::higher(self.nodes[position], old) {
                self.up(position);
            } else {
                self.down(position);
            }
        }
    }
    fn pop(&mut self) -> Option<(Reverse<u32>, usize)> {
        if self.nodes.is_empty() {
            None
        } else {
            let (degree, cell) = self.remove(0);
            Some((Reverse(degree), cell as usize))
        }
    }
    fn peek(&self) -> Option<(Reverse<u32>, usize)> {
        self.nodes
            .first()
            .map(|&(degree, cell)| (Reverse(degree), cell as usize))
    }
    fn update(&mut self, cell: usize, degree: u32) {
        if degree > 1 {
            self.push((Reverse(degree), cell));
        } else if self.location[cell] != u32::MAX {
            self.remove(self.location[cell] as usize);
        }
    }
}
