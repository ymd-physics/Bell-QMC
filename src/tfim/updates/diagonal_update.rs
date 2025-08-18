/*********************************************************************************
    For each spin, the state is
        (rz, rx) = (0,0), (0,1), (1,0), (1,1)
    -------------------------------------------------------------------
    The hash map considered here is 
        op_string[p] = 4 * b + t 
    then
        [-1, -1], op: -1, null operator
        [0, s], op = 4 * s + 0, diag, site operator 
        [1, s], op = 4 * s + 1, off-diag, site operator
        [2, b], op = 4 * b + 2, diag, bond operator
        [3, b], op = 4 * b + 3, off-diag, bond operator
    -------------------------------------------------------------------
    The legs are labelled as 
        2    3
       [][][][]
        0    1
    with "UP" to be the imaginary time direction
        The "left/right" qudits then save [0, 1] legs
    -------------------------------------------------------------------
    The updates on off-diag site and bond operators are individual
*********************************************************************************/
use crate::tfim::TFIModel;
use crate::aux::{NULL_OP, NULL_QUDIT};

impl TFIModel {
    pub fn diag_update(&mut self) {
        let mut op: i32;
        let mut remainder: usize;
        let mut new_bond: usize;
        let mut new_site: usize;
        let mut the_bond: usize;
        let mut the_site: usize;
        let mut the_qudit_left: u8;
        let mut the_qudit_right: u8;

        for p in 0..self.m {
            op = self.op_string[p];

            // ---------------------------------------
            //  Encounter a null operator
            // ---------------------------------------
            if op == NULL_OP {
                // First decide whether to insert
                let the_prob = self.add_factor / (self.m - self.n) as f64;
                if (the_prob >= 1.0) || (self.rand_prob() <= the_prob) {
                    // then decide which to insert
                    if self.rand_prob() < self.selection_prob {
                        new_site = self.rand_site();
                        self.op_string[p] = (4 * new_site) as i32;
                        self.n += 1;
    
                        the_qudit_left = self.qudits[new_site];
                        the_qudit_right = NULL_QUDIT;
                    }
                    
                    else {
                        new_bond = self.rand_bond();
                        self.op_string[p] = (4 * new_bond + 2) as i32;
                        self.n += 1;

                        the_qudit_left = self.qudits[self.b_sites[new_bond][0]];
                        the_qudit_right = self.qudits[self.b_sites[new_bond][1]];
                    }
                }

                else {
                    the_qudit_left = NULL_QUDIT;
                    the_qudit_right = NULL_QUDIT;
                }
            }

            else {
                remainder = (op % 4) as usize;

                match remainder {
                    // ---------------------------------------
                    //  Encounter a diagonal site operator
                    // ---------------------------------------
                    0 => {
                        let the_prob = self.remove_factor * (self.m - self.n + 1) as f64;
                        if (the_prob >= 1.0) || (self.rand_prob() <= the_prob)  {
                            self.op_string[p] = NULL_OP;
                            self.n -= 1;

                            the_qudit_left = NULL_QUDIT;
                            the_qudit_right = NULL_QUDIT;
                        }

                        else {
                            the_site = (op / 4) as usize;
                            the_qudit_left = self.qudits[the_site];
                            the_qudit_right = NULL_QUDIT;
                        }
                    }

                    // ---------------------------------------
                    //  Encounter a diagonal bond operator
                    // ---------------------------------------
                    2 => {
                        let the_prob = self.remove_factor * (self.m - self.n + 1) as f64;
                        if (the_prob >= 1.0) || (self.rand_prob() <= the_prob)  {
                            self.op_string[p] = NULL_OP;
                            self.n -= 1;
                            the_qudit_left = NULL_QUDIT;
                            the_qudit_right = NULL_QUDIT;
                        }

                        else {
                            the_bond = (op / 4) as usize;
                            the_qudit_left = self.qudits[self.b_sites[the_bond][0]];
                            the_qudit_right = self.qudits[self.b_sites[the_bond][1]];
                        }
                    }

                    // --------------------------------------------
                    //  Encounter an off-diagonal site operator
                    // --------------------------------------------
                    1 => {
                        the_site = (op / 4) as usize;
                        the_qudit_left = self.qudits[the_site];
                        the_qudit_right = NULL_QUDIT;

                        // The site off-diagonal operator let
                        //      rz, rx --> rz, (rx + 1) mod 2
                        self.qudits[the_site] ^= 0b01;
                    }

                    // --------------------------------------------
                    //  Encounter an off-diagonal bond operator
                    // --------------------------------------------
                    3 => {
                        the_bond = (op / 4) as usize;
                        the_qudit_left = self.qudits[self.b_sites[the_bond][0]];
                        the_qudit_right = self.qudits[self.b_sites[the_bond][1]];

                        // The bond off-diagonal operator let
                        //      rz, rx --> (rz + 1) mod 2, rx 
                        // for the two sites it acts on
                        self.qudits[self.b_sites[the_bond][0]] ^= 0b10;
                        self.qudits[self.b_sites[the_bond][1]] ^= 0b10;
                    }

                    // --------------------------------------------
                    //  Others (errors)
                    // --------------------------------------------
                    _ => {
                        panic!("Invalid bond occur with op_string[p] = {}", p);
                    }
                }
            }

            self.left_qudits[p] = the_qudit_left;
            self.right_qudits[p] = the_qudit_right;
        }
    }

    pub fn refresh_left_right_qudits(&mut self) {
        let mut op: i32;
        let mut remainder: usize;

        let mut the_bond: usize;
        let mut the_site: usize;
        let mut the_qudit_left: u8;
        let mut the_qudit_right: u8;

        for p in 0..self.m {
            op = self.op_string[p];

            // ---------------------------------------
            //  Encounter a null operator
            // ---------------------------------------
            if op == NULL_OP {
                the_qudit_left = NULL_QUDIT;
                the_qudit_right = NULL_QUDIT;    
            }

            else {
                remainder = (op % 4) as usize;

                match remainder {
                    // ---------------------------------------
                    //  Encounter a diagonal site operator
                    // ---------------------------------------
                    0 => {
                        the_site = (op / 4) as usize;
                        the_qudit_left = self.qudits[the_site];
                        the_qudit_right = NULL_QUDIT;
                    }

                    // ---------------------------------------
                    //  Encounter a diagonal bond operator
                    // ---------------------------------------
                    2 => {
                        the_bond = (op / 4) as usize;
                        the_qudit_left = self.qudits[self.b_sites[the_bond][0]];
                        the_qudit_right = self.qudits[self.b_sites[the_bond][1]];
                    }

                    // --------------------------------------------
                    //  Encounter an off-diagonal site operator
                    // --------------------------------------------
                    1 => {
                        the_site = (op / 4) as usize;
                        the_qudit_left = self.qudits[the_site];
                        the_qudit_right = NULL_QUDIT;

                        // The site off-diagonal operator let
                        //      rz, rx --> rz, (rx + 1) mod 2
                        self.qudits[the_site] ^= 0b01;
                    }

                    // --------------------------------------------
                    //  Encounter an off-diagonal bond operator
                    // --------------------------------------------
                    3 => {
                        the_bond = (op / 4) as usize;
                        the_qudit_left = self.qudits[self.b_sites[the_bond][0]];
                        the_qudit_right = self.qudits[self.b_sites[the_bond][1]];

                        // The bond off-diagonal operator let
                        //      rz, rx --> (rz + 1) mod 2, rx 
                        // for the two sites it acts on
                        self.qudits[self.b_sites[the_bond][0]] ^= 0b10;
                        self.qudits[self.b_sites[the_bond][1]] ^= 0b10;

                    }

                    // --------------------------------------------
                    //  Others (errors)
                    // --------------------------------------------
                    _ => {
                        panic!("Invalid bond occur with op_string[p] = {}", p);
                    }
                }
            }

            self.left_qudits[p] = the_qudit_left;
            self.right_qudits[p] = the_qudit_right;
        }
    }

    pub fn diag_update_with_measure(&mut self) {
        let mut op: i32;
        let mut remainder: usize;
        let mut new_bond: usize;
        let mut new_site: usize;
        let mut the_bond: usize;
        let mut the_site: usize;
        let mut the_qudit_left: u8;
        let mut the_qudit_right: u8;

        for p in 0..self.m {
            // ---------------------------------------
            //  Measurements on this time slice
            // ---------------------------------------
            self.measure();

            // ---------------------------------------
            //  Encounter a null operator
            // ---------------------------------------
            op = self.op_string[p];

            if op == NULL_OP {
                // First decide whether to insert
                let the_prob = self.add_factor / (self.m - self.n) as f64;
                if (the_prob >= 1.0) || (self.rand_prob() <= the_prob) {
                    // then decide which to insert
                    if self.rand_prob() < self.selection_prob {
                        new_site = self.rand_site();
                        self.op_string[p] = (4 * new_site) as i32;
                        self.n += 1;
    
                        the_qudit_left = self.qudits[new_site];
                        the_qudit_right = NULL_QUDIT;
                    }
                    
                    else {
                        new_bond = self.rand_bond();
                        self.op_string[p] = (4 * new_bond + 2) as i32;
                        self.n += 1;

                        the_qudit_left = self.qudits[self.b_sites[new_bond][0]];
                        the_qudit_right = self.qudits[self.b_sites[new_bond][1]];
                    }
                }

                else {
                    the_qudit_left = NULL_QUDIT;
                    the_qudit_right = NULL_QUDIT;
                }
            }

            else {
                remainder = (op % 4) as usize;

                match remainder {
                    // ---------------------------------------
                    //  Encounter a diagonal site operator
                    // ---------------------------------------
                    0 => {
                        let the_prob = self.remove_factor * (self.m - self.n + 1) as f64;
                        if (the_prob >= 1.0) || (self.rand_prob() <= the_prob)  {
                            self.op_string[p] = NULL_OP;
                            self.n -= 1;

                            the_qudit_left = NULL_QUDIT;
                            the_qudit_right = NULL_QUDIT;
                        }

                        else {
                            the_site = (op / 4) as usize;
                            the_qudit_left = self.qudits[the_site];
                            the_qudit_right = NULL_QUDIT;
                        }
                    }

                    // ---------------------------------------
                    //  Encounter a diagonal bond operator
                    // ---------------------------------------
                    2 => {
                        let the_prob = self.remove_factor * (self.m - self.n + 1) as f64;
                        if (the_prob >= 1.0) || (self.rand_prob() <= the_prob)  {
                            self.op_string[p] = NULL_OP;
                            self.n -= 1;
                            the_qudit_left = NULL_QUDIT;
                            the_qudit_right = NULL_QUDIT;
                        }

                        else {
                            the_bond = (op / 4) as usize;
                            the_qudit_left = self.qudits[self.b_sites[the_bond][0]];
                            the_qudit_right = self.qudits[self.b_sites[the_bond][1]];
                        }
                    }

                    // --------------------------------------------
                    //  Encounter an off-diagonal site operator
                    // --------------------------------------------
                    1 => {
                        the_site = (op / 4) as usize;
                        the_qudit_left = self.qudits[the_site];
                        the_qudit_right = NULL_QUDIT;

                        // The site off-diagonal operator let
                        //      rz, rx --> rz, (rx + 1) mod 2
                        self.qudits[the_site] ^= 0b01;
                    }

                    // --------------------------------------------
                    //  Encounter an off-diagonal bond operator
                    // --------------------------------------------
                    3 => {
                        the_bond = (op / 4) as usize;
                        the_qudit_left = self.qudits[self.b_sites[the_bond][0]];
                        the_qudit_right = self.qudits[self.b_sites[the_bond][1]];

                        // The bond off-diagonal operator let
                        //      rz, rx --> (rz + 1) mod 2, rx 
                        // for the two sites it acts on
                        self.qudits[self.b_sites[the_bond][0]] ^= 0b10;
                        self.qudits[self.b_sites[the_bond][1]] ^= 0b10;
                    }

                    // --------------------------------------------
                    //  Others (errors)
                    // --------------------------------------------
                    _ => {
                        panic!("Invalid bond occur with op_string[p] = {}", p);
                    }
                }
            }

            self.left_qudits[p] = the_qudit_left;
            self.right_qudits[p] = the_qudit_right;
        }
    }
}