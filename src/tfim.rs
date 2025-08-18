use prng_mt::mt19937::MT19937;
pub mod init;
pub mod random;
pub mod updates;
pub mod measure;
pub mod stack;

pub struct TFIModel {
    // ----------------------------------------------------------------
    //  Basic params
    // ----------------------------------------------------------------
    l: usize,
    beta: f64,

    j: f64,     // strength of ZZ couplings
    h: f64,     // strength of the external fields

    n: usize,       // number of null operators
    pub m: usize,       // truncation order of the series

    // --------------------------------------------------------
    //  Lattice
    // --------------------------------------------------------
    pub num_sites: usize,             // number of sites
    num_bonds: usize,                 // number of bonds
    b_sites: Vec<Vec<usize>>,         // record the linking of the lattice

    // --------------------------------------------------------
    //  Some frequently-used factors
    // --------------------------------------------------------
    selection_prob: f64,
    add_factor: f64,
    remove_factor: f64,

    // -------------------------------------------------
    //  PRNG with MT19937
    // -------------------------------------------------
    rng: MT19937,

    // -----------------------------------------------
    //  Data structures for configuration updates
    // -----------------------------------------------
    qudits: Vec<u8>,         // the qudits at the zero time
    left_qudits: Vec<u8>,    // the qudits for the operator at each time slice
    right_qudits: Vec<u8>,   
    op_string: Vec<i32>,
    v_first: Vec<i32>,
    v_last: Vec<i32>,
    vertex_list: Vec<i32>,

    // --------------------------------------------------------------------
    //  Two internal stacks: Pre-allocating memory makes it faster
    // --------------------------------------------------------------------
    stack: Vec<usize>,
    top: i32,
    flip: i32,

    // ----------------------------------------
    //  Measurements related
    // ----------------------------------------
    subsystem: Vec<usize>,
    system: Vec<usize>,
    pub purity: f64,
    pub partial_purity: f64, 
    pub zz_corr_2: Vec<f64>,    // record the correlations between the 0th site and others
    pub xx_corr_2: Vec<f64>, 

    // ----------------------------------------
    //  For testing the parity
    // ----------------------------------------
    pub num_parity_odd: usize,
    pub num_parity_even: usize,
}