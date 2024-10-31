use std::env;

use log_stats::logstats::simple;

fn main() {
    
    let args: Vec<String> = env::args().collect();
    //arg[0] is the output package directory
    //arg[1] should be the log file path.
    if args.len() == 1 {
        print!("missing the log file path\n");
    } else if args.len() != 2 {
        print!("expecting only 1 argument, the log file path\n");
    } else {
        print!("reading file {}", args[0]);
        simple(args[1].clone());
    }
}
