use prng_mt::MT19937;
use crate::tfim::TFIModel;
use crate::aux::{NULL_OP, NULL_QUDIT, EMPTY};

impl TFIModel {
    pub fn new(para_l: usize, para_beta: f64, para_j: f64, para_h: f64, para_seed: u32) -> Self {
        Self {
            // ----------------------------------------------------------------
            //  Basic params
            // ----------------------------------------------------------------
            l: para_l,
            beta: para_beta,
            j: para_j,
            h: para_h,
            n: 0,
            m: 10,

            // --------------------------------------------------------
            //  Lattice
            // --------------------------------------------------------
            num_sites: para_l,
            num_bonds: 0,                 // default
            b_sites: Vec::new(),         // default

            // --------------------------------------------------------
            //  Some frequently-used factors
            // --------------------------------------------------------
            selection_prob: 0.0,    // default
            add_factor: 0.0,        // default
            remove_factor: 0.0,     // default

            // ----------------------------------------
            //  Random number generator
            // ----------------------------------------
            rng: MT19937::new(para_seed),

            // ----------------------------------------
            //  Data structures for configuration
            // ----------------------------------------
            qudits: Vec::new(),
            left_qudits: Vec::new(),
            right_qudits: Vec::new(),

            op_string: Vec::new(),
            v_first: Vec::new(),
            v_last: Vec::new(),
            vertex_list: Vec::new(),
    
            // ------------------------------------
            //  Internal stacks
            // ------------------------------------
            stack: Vec::new(),
            top: 0,
            flip: 0,

            // ----------------------------------------
            //  Measurements related
            // ----------------------------------------
            system: Vec::new(),
            subsystem: Vec::new(),
            partial_purity: 0.0, 
            purity: 0.0,
            zz_corr_2: vec![0.0; para_l],
            xx_corr_2: vec![0.0; para_l],

            // ----------------------------------------
            //  For testing the parity
            // ----------------------------------------
            num_parity_odd: 0,
            num_parity_even: 0,
        }
    }

    pub fn init(&mut self) {
        self.num_bonds = self.l - 1;

        // ------------------------------------
        //  Make the lattice (OBC)
        // ------------------------------------


        self.b_sites = {
            let mut lattice: Vec<Vec<usize>> = Vec::new();
            //  make the 1D lattice (OBC)
            for i in 0..(self.num_sites - 1) {
                lattice.push(vec![i, i + 1]);
            }
            lattice.clone()
        };

        // --------------------------------------------
        //  Initialize the frequently-used factors
        // --------------------------------------------
        self.selection_prob = self.h * self.num_sites as f64 / (
            self.h * self.num_sites as f64 + self.j * self.num_bonds as f64 
        );

        self.add_factor = self.beta * (
            self.h * self.num_sites as f64 + self.j * self.num_bonds as f64
        );
        self.remove_factor = 1.0 / self.add_factor;

        // --------------------------------------------
        //  Initialize the initial states
        // --------------------------------------------
        // self.qudits = (0..self.num_sites)
        //     .map(|_| self.rand_qudit())
        //     .collect();
        self.qudits = vec![0; self.num_sites];

         // --------------------------------------------------------------------
        //  Initialize the data structures for saving states at time "p"
        // --------------------------------------------------------------------
        self.left_qudits = vec![NULL_QUDIT; self.m];
        self.right_qudits = vec![NULL_QUDIT; self.m];

        // --------------------------------------------------------
        //  Initialize data structures related to operator string
        // --------------------------------------------------------
        self.op_string = vec![NULL_OP; self.m];    

        self.v_first = vec![EMPTY; self.num_sites];     // for considering the two virtual bonds 
        self.v_last = vec![EMPTY; self.num_sites];
        self.vertex_list = vec![EMPTY; 4 * self.m];
  
        // ------------------------------------------------------------------------
        // Initialize for the internal stack (the capacity should be large enough)
        // ------------------------------------------------------------------------
        self.stack = vec![0; 8 * self.m];

        // -------------------------------------------
        //  Decide the environment to be traced out
        // -------------------------------------------
        self.system = (0..self.num_sites).collect();
        self.subsystem = (0..(self.num_sites / 2)).collect();
    }
}