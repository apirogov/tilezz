pub struct SparseTable {
    table: Vec<Vec<usize>>,
    log: Vec<usize>,
}

impl SparseTable {
    pub fn new(data: &[usize]) -> Self {
        let n = data.len();
        if n == 0 {
            return SparseTable {
                table: vec![],
                log: vec![],
            };
        }

        let max_k = floor_log2(n);
        let mut table = Vec::with_capacity(max_k + 1);

        table.push(data.to_vec());

        for k in 1..=max_k {
            let half = 1 << (k - 1);
            let prev = &table[k - 1];
            let len = n.saturating_sub(half);
            let mut row = vec![0usize; n];
            for i in 0..len {
                row[i] = prev[i].min(prev[i + half]);
            }
            row[len..n].copy_from_slice(&prev[len..n]);
            table.push(row);
        }

        let mut log = vec![0usize; n + 1];
        for i in 2..=n {
            log[i] = log[i / 2] + 1;
        }

        SparseTable { table, log }
    }

    pub fn query(&self, l: usize, r: usize) -> usize {
        assert!(l <= r, "query requires l <= r");
        let k = self.log[r - l + 1];
        self.table[k][l].min(self.table[k][r + 1 - (1 << k)])
    }
}

fn floor_log2(n: usize) -> usize {
    if n == 0 {
        return 0;
    }
    (usize::BITS - 1 - n.leading_zeros()) as usize
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_rmq(data: &[usize], l: usize, r: usize) -> usize {
        data[l..=r].iter().copied().min().unwrap()
    }

    #[test]
    fn empty() {
        let st = SparseTable::new(&[]);
        assert!(st.table.is_empty());
    }

    #[test]
    fn single_element() {
        let st = SparseTable::new(&[42]);
        assert_eq!(st.query(0, 0), 42);
    }

    #[test]
    fn all_same() {
        let data = vec![5, 5, 5, 5, 5];
        let st = SparseTable::new(&data);
        for l in 0..data.len() {
            for r in l..data.len() {
                assert_eq!(st.query(l, r), 5);
            }
        }
    }

    #[test]
    fn ascending() {
        let data: Vec<usize> = (0..10).collect();
        let st = SparseTable::new(&data);
        for l in 0..data.len() {
            for r in l..data.len() {
                assert_eq!(st.query(l, r), l, "rmq({l}, {r})");
            }
        }
    }

    #[test]
    fn descending() {
        let data: Vec<usize> = (0..10).rev().collect();
        let st = SparseTable::new(&data);
        for l in 0..data.len() {
            for r in l..data.len() {
                assert_eq!(st.query(l, r), naive_rmq(&data, l, r), "rmq({l}, {r})");
            }
        }
    }

    #[test]
    fn all_pairs_vs_naive() {
        let data = vec![3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5];
        let st = SparseTable::new(&data);
        for l in 0..data.len() {
            for r in l..data.len() {
                assert_eq!(st.query(l, r), naive_rmq(&data, l, r), "rmq({l}, {r})");
            }
        }
    }

    #[test]
    fn power_of_two_length() {
        let data: Vec<usize> = (0..16).collect();
        let st = SparseTable::new(&data);
        assert_eq!(st.query(0, 15), 0);
        assert_eq!(st.query(7, 8), 7);
        assert_eq!(st.query(15, 15), 15);
    }

    #[test]
    fn random_vs_naive() {
        let mut s: u64 = 31415;
        let data: Vec<usize> = (0..100)
            .map(|_| {
                s = s
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                (s % 1000) as usize
            })
            .collect();
        let st = SparseTable::new(&data);

        let mut s2: u64 = 27182;
        for _ in 0..500 {
            s2 = s2
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let l = (s2 as usize) % data.len();
            s2 = s2
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let r = l + ((s2 as usize) % (data.len() - l));
            assert_eq!(st.query(l, r), naive_rmq(&data, l, r), "rmq({l}, {r})");
        }
    }
}
