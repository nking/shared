use std::f64;
// cannot declare struct fields modifiable,
// but can modify them with extended function mutable self ref
pub struct StatsHolder {
    x_sum: u64,
    x_sum_sq: u64,
    n : u32,
    _description: String,
}

// implement the behavior
impl StatsHolder {
    // Not get_first.
    pub fn accumulate(& mut self, x : u64) -> () {
        self.x_sum += x;
        self.x_sum_sq += x*x;
        self.n += 1;
    }

    // Not get_first_mut, get_mut_first, or mut_first.
    pub fn calc_mean_stdev(& self) -> (f64, f64) {
        let m : f64 = (self.x_sum as f64) / (self.n as f64);
        let stdev : f64 = (((self.x_sum_sq as f64) - m * m * (self.n as f64))/(self.n as f64 - 1.0f64)).sqrt();

        (m, stdev)
    }

    pub fn is_empty(& self) -> bool {
        return self.n == 0;
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
