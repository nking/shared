use std::{f64, u64};
// cannot declare struct fields modifiable,
// but can modify them with extended function mutable self ref
pub struct StatsHolder {
    x_sum: u64,
    x_sum_sq: u64,
    n : u32,
    _description: String,
}

pub struct StatsHolderMed {
    x: Vec<u64>,
    stats : StatsHolder,
}

// implement the behavior
impl StatsHolder {
    pub fn accumulate(& mut self, x : u64) -> () {
        self.x_sum += x;
        self.x_sum_sq += x*x;
        self.n += 1;
    }

    pub fn calc_mean_stdev(& self) -> (f64, f64) {
        let m : f64 = (self.x_sum as f64) / (self.n as f64);
        let stdev : f64 = (((self.x_sum_sq as f64) - m * m * (self.n as f64))/(self.n as f64 - 1.0f64)).sqrt();

        (m, stdev)
    }

    pub fn is_empty(& self) -> bool {
        return self.n == 0;
    }
}

impl StatsHolderMed {
    pub fn accumulate(& mut self, x : u64) -> () {
        self.stats.accumulate(x);
        self.x.push(x);
    }

    pub fn calc_mean_stdev(& self) -> (f64, f64) {
        self.stats.calc_mean_stdev()
    }

    pub fn calc_median(& mut self) -> f64 {
        println!("{:?} x len={:?}", self.stats._description, self.x.len());
        //TODO: consider making a crate for the linear time median algorithms.
        //see src/main/java/algorithms/MedianOfMediansSelect
        if self.x.is_empty() {
            return 0.0f64;
        }
        if !self.x.is_sorted() {
            self.x.sort();
        }
        if (self.x.len() & 1) != 1 {
            return self.x[self.x.len()/2].clone() as f64;
        } else {
            let m1 = self.x[self.x.len()/2];
            let m2 = self.x[(self.x.len()/2) - 1];
            return (m1 as f64 + m2 as f64)/2.0;
        }
    }

    pub fn is_empty(& self) -> bool {
        return self.x.is_empty();
    }
}

pub fn build_stats_holder(description: String) -> StatsHolder {
    StatsHolder {
        _description: description,
        x_sum: 0,
        x_sum_sq: 0,
        n: 0,
    }
}

pub fn build_stats_holder_med(description: String) -> StatsHolderMed {
    let stats = StatsHolder {
        _description: description,
        x_sum: 0,
        x_sum_sq: 0,
        n: 0,
    };

    StatsHolderMed {
        x: Vec::<u64>::new(),
        stats: stats,
    }
}
