use crate::tfim::TFIModel;
use crate::aux::{NULL_OP, EMPTY};


impl TFIModel {
    pub fn make_vertex_list(&mut self) {
        let mut op: i32;
        let mut b_p: usize;
        let mut v_leg0: i32;
        let mut s0: usize;
        let mut s1: usize;
        let mut s0_v_last: i32;
        let mut s1_v_last: i32;

        // ====================================================
        //  Initialize "vertex_list", "v_first", "v_last"
        // ====================================================
        for v in 0..4 * self.m {
            self.vertex_list[v] = EMPTY;
        }
            
        for s in 0..self.num_sites {
            self.v_first[s] = EMPTY;
            self.v_last[s] = EMPTY;
        }

        // ===================================
        //  Make the vertex_list
        // ===================================
        for p in 0..self.m {
            op = self.op_string[p];
            
            if op != NULL_OP {
                // -.-.-.-.-.-.-.-.-.
                //  Bond operator
                // -.-.-.-.-.-.-.-.-.
                if op % 4 >= 2 {
                    b_p = (op / 4) as usize;
                    s0 = self.b_sites[b_p][0];
                    s1 = self.b_sites[b_p][1];
                    v_leg0 = (4 * p) as i32;

                    s0_v_last = self.v_last[s0];
                    s1_v_last = self.v_last[s1];

                    if s0_v_last > EMPTY {
                        self.vertex_list[s0_v_last as usize] = v_leg0;
                        self.vertex_list[v_leg0 as usize] = s0_v_last;
                    } 
                    
                    else {
                        self.v_first[s0] = v_leg0;
                    }
                        
                    self.v_last[s0] = v_leg0 + 2;

                    if s1_v_last > EMPTY {
                        self.vertex_list[s1_v_last as usize] = v_leg0 + 1;
                        self.vertex_list[(v_leg0 + 1) as usize] = s1_v_last;
                    } 
                    
                    else {
                        self.v_first[s1] = v_leg0 + 1;
                    }
                        
                    self.v_last[s1] = v_leg0 + 3;
                }

                // -.-.-.-.-.-.-.-.-.
                //  Site operator
                // -.-.-.-.-.-.-.-.-.
                else {
                    s0 = (op / 4) as usize;
                    v_leg0 = (4 * p) as i32;
                    s0_v_last = self.v_last[s0];

                    if s0_v_last > EMPTY {
                        self.vertex_list[s0_v_last as usize] = v_leg0;
                        self.vertex_list[v_leg0 as usize] = s0_v_last;
                    } 

                    else {
                        self.v_first[s0] = v_leg0;
                    }
                        
                    self.v_last[s0] = v_leg0 + 2;
                }
            }
        }

        // ===================================
        //  Make PBC correction
        // ===================================
        let mut s_v_first;
        let mut s_v_last;

        for s in 0..self.num_sites {
            s_v_first = self.v_first[s];
            
            if s_v_first != EMPTY {
                s_v_last = self.v_last[s];
                self.vertex_list[s_v_first as usize] = s_v_last;
                self.vertex_list[s_v_last as usize] = s_v_first;
            }
        }
    }

    pub fn make_dual_vertex_list(&mut self) {
        // Only the off-diagonal site operator is stretched
        let mut op: i32;
        let mut b_p: usize;
        let mut remainder: usize;

        let mut v_leg0: i32;
        let mut b0: usize;
        let mut b1: usize;

        let mut b0_v_last: i32;
        let mut b1_v_last: i32;

        // ====================================================
        //  Initialize "vertex_list", "v_first", "v_last"
        // ====================================================
        for v in 0..4 * self.m {
            self.vertex_list[v] = EMPTY;
        }
            
        for s in 0..self.num_sites {
            self.v_first[s] = EMPTY;
            self.v_last[s] = EMPTY;
        }

        // ===================================
        //  Make the vertex_list
        // ===================================
        for p in 0..self.m {
            op = self.op_string[p];

            if op != NULL_OP {
                b_p = (op / 4) as usize;
                remainder = (op % 4) as usize;
                v_leg0 = (4 * p) as i32;

                match remainder {
                    2 | 3 => {
                        // -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                        //  Bond operator
                        // -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                        b0 = self.b_sites[b_p][0];  // use the label of the left site for labelling the bond
                        b0_v_last = self.v_last[b0];

                        if b0_v_last > EMPTY {
                            self.vertex_list[b0_v_last as usize] = v_leg0;
                            self.vertex_list[v_leg0 as usize] = b0_v_last;
                        } 

                        else {
                            self.v_first[b0] = v_leg0;
                        }

                        self.v_last[b0] = v_leg0 + 1;
                    }

                    0 | 1 => {
                        // -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                        //  Off-diagonal site operator
                        // -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
                        b0 = (b_p - 1 + self.num_sites) % self.num_sites;
                        b1 = b_p;

                        {
                            b0_v_last = self.v_last[b0 as usize];
                            if b0_v_last > EMPTY {
                                self.vertex_list[b0_v_last as usize] = v_leg0;
                                self.vertex_list[v_leg0 as usize] = b0_v_last;
                            } 

                            else {
                                self.v_first[b0 as usize] = v_leg0;
                            }

                            self.v_last[b0 as usize] = v_leg0 + 1;
                        }

                        {
                            b1_v_last = self.v_last[b1];
                            if b1_v_last > EMPTY {
                                self.vertex_list[b1_v_last as usize] = v_leg0 + 2;
                                self.vertex_list[(v_leg0 + 2) as usize] = b1_v_last;
                            } 

                            else {
                                self.v_first[b1] = v_leg0 + 2;
                            }

                            self.v_last[b1] = v_leg0 + 3;
                        }
                    }
                    _ => { panic!(); }
                }
            }
        }

        // -===================================
        //  Make PBC correction
        // ===================================
        let mut s_v_first;
        let mut s_v_last;

        for s in 0..self.num_sites {
            s_v_first = self.v_first[s];
            
            if s_v_first != EMPTY {
                s_v_last = self.v_last[s];
                self.vertex_list[s_v_first as usize] = s_v_last;
                self.vertex_list[s_v_last as usize] = s_v_first;
            }
        }
    }

}