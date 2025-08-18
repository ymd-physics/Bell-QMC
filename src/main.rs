/*************************************************************************************
 *  SSE for 1D TFIM OBC under Bell basis, GS version
 *  Author: Yi-Ming Ding
 *  Updated: Mar 12, 2025
 ************************************************************************************/
pub mod aux;
pub mod tfim; 
use std::{env, fs::File, io::Write, time::Instant};

fn main() {
    let start_time = Instant::now();

    // ========================================================
    //  Collect params from shell
    // ========================================================
    let args: Vec<String> = env::args().collect();
    let para_l: usize = args[1].parse().unwrap();
    let mut para_beta: f64 = args[2].parse().unwrap();
    let para_j: f64 = args[3].parse().unwrap();
    let para_h: f64 = args[4].parse().unwrap();
    let num_thm: usize = args[5].parse().unwrap();
    let num_stat: usize = args[6].parse().unwrap();
    let num_bins: usize = args[7].parse().unwrap();
    let target_dir: String = args[8].parse().unwrap();
    let para_seed: u32 = args[9].parse().unwrap();

    // ===============================================================
    //  Report the environment
    // ================================================================
    aux::print_horizontal_line(77, "-");
    println!("■ Bell-QMC for 1D TFIM (OBC, ground state simulation)");
    println!("■ l = {para_l}, beta = {para_beta}, J = {para_j}, h = {para_h}");
    println!("■ num_thm = {num_thm}, num_stat = {num_stat}, num_bins = {num_bins}, seed = {para_seed}");
    para_beta *= 2.0;

    // ===============================================================
    //  Preparing for writing the results
    // ===============================================================
    let mut file_purity: File = aux::create_new_file(format!("{}/purity.dat", target_dir)); 
    let mut file_renyi_ee: File = aux::create_new_file(format!("{}/renyi2_ee.dat", target_dir)); 
    let mut file_zz: File = aux::create_new_file(format!("{}/zz_corr_2.dat", target_dir));
    let mut file_xx: File = aux::create_new_file(format!("{}/xx_corr_2.dat", target_dir));

    // ===============================================================
    //  Monte Carlo simulations
    // ===============================================================
    let mut model = tfim::TFIModel::new(para_l,para_beta, para_j, para_h, para_seed);
    model.init();

    println!("\t---> Thermalizing...");
    for _ in 0..num_thm { 
        model.mc_thermalizing(); 
    }

    // notice that we utilize samples in time slices, thus "num_stat" can be modified
    println!("\t---> Maximum cut-off = {}", model.m);
    let num_samples: f64 = num_stat as f64 * model.m as f64;    // use f64 to avoid the overflow
    println!("\t---> Total number of samples = {} * {} = {} ", num_stat, model.m, num_samples);
    println!("\t---> Sampling and measuring...");
    for b in 0..num_bins {
        println!("\t\t# bin {}...", b);
        model.ini_measure();
        for _ in 0..num_stat {
            model.mc_sampling();
        }
        model.statisticize(num_samples); 

        // ------------------------------------
        //  Saving the data
        // ------------------------------------
        file_purity.write_all(format!("{:<16.10}\n", model.purity).as_bytes()).unwrap();
        file_renyi_ee.write_all(format!("{:<16.10}\n", -1.0 * model.partial_purity.ln()).as_bytes()).unwrap();
        file_zz.write_all(format!("{}\n",
                                  model.zz_corr_2.iter().map(|x| format!("{:<16.10}", x)).collect::<Vec<_>>().join("\t")
        ).as_bytes()).unwrap();
        file_xx.write_all(format!("{}\n",
                                  model.xx_corr_2.iter().map(|x| format!("{:<16.10}", x)).collect::<Vec<_>>().join("\t")
        ).as_bytes()).unwrap();
    }
    
    // =============================================
    //  Report the runtime
    // =============================================
    aux::print_horizontal_line(77, "-");
    aux::report_time(start_time);
}