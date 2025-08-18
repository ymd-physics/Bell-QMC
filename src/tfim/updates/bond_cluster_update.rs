use crate::tfim::TFIModel;
use crate::{flip_operator, flip_rz, flip_rx, get_rx};
use crate::aux::{EMPTY, FLIPPED, FREE_SPIN, NOT_FLIPPED};

#[inline]
fn to_back(v: usize) -> usize { v ^ 01 }

impl TFIModel {
    fn link_to_valid_dual_cluster_leg(&mut self, v: usize) -> i32 {
        // --------------------------------------------------
        //  It stops at 
        //      (i) a flippable bond operator
        //      (ii) an off-diagonal site operator
        // --------------------------------------------------
        if self.vertex_list[v] < 0 {
            self.vertex_list[v]
        }

        else {
            let mut v0 : usize = v;     // this is the start point 
            let mut v1: usize;
            let mut the_p1: usize;
            let mut op1: i32; 
            let mut remainder1: usize;

            loop {
                v1 = self.vertex_list[v0] as usize;
                the_p1 = v1 / 4;
                op1 = self.op_string[the_p1];
                remainder1 = (op1 % 4) as usize;
                
                // an off-diagonal site is always okay
                // A diag-bond operator that dose not satisfy the constraint is not valid
                if (remainder1 == 0) || (remainder1 == 2 && get_rx!(self.left_qudits[the_p1]) != get_rx!(self.right_qudits[the_p1])) {
                    v0 = to_back(v1); 
                }
                else {
                    break;
                }
            }

            v1 as i32 
        }
    }

    pub fn bond_cluster_update(&mut self) {
        self.make_dual_vertex_list();

        let mut op: i32;
        let mut the_p: usize;
        let mut remainder: usize;

        let mut v0: usize;
        self.stack_initialize();        // for growing the cluster

        for v in 0..4 * self.m {
            if self.vertex_list[v] < 0 { 
                continue; 
            }

            // ------------------------------------------------------------------------------
            // We must ensure it is a valid leg, i.e. it belongs to a valid bond operator 
            // ------------------------------------------------------------------------------
            the_p = v / 4;
            op = self.op_string[the_p];
            remainder = (op % 4) as usize;

            // --------------------------------------------------------------------
            // an off-diagonal bond operator is always flippable
            // But for a diag-bond operator, we skip it if r^x != r^x
            // --------------------------------------------------------------------
            if remainder >= 2 {
                if (remainder == 3) || (get_rx!(self.left_qudits[the_p]) == get_rx!(self.right_qudits[the_p])) {
                    self.flip = if self.rand_prob() > 0.5 { FLIPPED } else { NOT_FLIPPED };
                    self.stack_push(v);
                    
                    // --------------------------------------------
                    //  Check whether "v" acorss no site op
                    // --------------------------------------------
                    let v1: usize = self.link_to_valid_dual_cluster_leg(v) as usize;
                    let the_p1: usize = v1 / 4;
                    let op1: i32 = self.op_string[the_p1];
                    let remainder1 = (op1 % 4) as usize;

                    if remainder1 >= 2 {
                        if self.flip == FLIPPED {
                            self.op_string[the_p] = flip_operator!(self.op_string[the_p]);
                            self.op_string[the_p1] = flip_operator!(self.op_string[the_p1]);
                        }
                        self.vertex_list[v] = self.flip;
                        self.vertex_list[v1] = self.flip;
                    }

                    // --------------------------------------------
                    //  Grow the dual cluster starting with "v"
                    // --------------------------------------------
                    else {
                        loop {
                            if self.top == EMPTY { break; }
                            self.make_bond_cluster();
                        }
                    }
                }
            }
        }

        // ============================================================================
        //  Update the qudits at zero time after the updates of bond operators
        // ============================================================================
        for b in 0..self.num_sites {
            if self.v_first[b] != FREE_SPIN {
                let v = self.v_first[b] as usize;  
                v0 = v;     // this is the start point

                loop {
                    if self.vertex_list[v0] < 0 {
                        if self.vertex_list[v0] == FLIPPED {
                            let s_left = b;
                            let s_right = (b + 1) % self.num_sites;
                            self.qudits[s_left] = flip_rz!(self.qudits[s_left]);
                            self.qudits[s_right] = flip_rz!(self.qudits[s_right]);
                        }
                        break;
                         
                    }
                    
                    else {
                        v0 = self.vertex_list[to_back(v0)] as usize;    
                        if v0 == v {
                            break;
                        }
                    }
                }
            }

            else {
                // Randomly change those isolated bonds (without changing the sector) 
                if self.rand_prob() > 0.5 {
                    let s_left = b;
                    let s_right = (b + 1) % self.num_sites;
                    self.qudits[s_left] = flip_rz!(self.qudits[s_left]);
                    self.qudits[s_right] = flip_rz!(self.qudits[s_right]);
                }
            }
        }

        // ============================================================================
        //  Flipping all the r^x globally
        // ============================================================================
        if self.rand_prob() > 0.5 {
            for s in 0..self.num_sites {
                self.qudits[s] = flip_rx!(self.qudits[s]);
            }
        }
    }

    fn make_bond_cluster(&mut self) {
        let the_p: usize;
        let op: i32;
        let remainder: usize;

        let v_start: usize = self.stack_pop();
        let v1: i32 = self.link_to_valid_dual_cluster_leg(v_start);

        if self.vertex_list[v_start] < 0 {
            /* ------------------------------------------------------------------------------------
                This means "vs" has been PROCESSED, so we do not need to consider it
                so then no need to process "vs" or link it to other vertices, just return.
            ------------------------------------------------------------------------------------ */
            return;
        }

        else if self.vertex_list[v1 as usize] != self.flip {
            /* ---------------------------------------------------------------------------------------
                This means "v1" has not been processed yet here, and we can push "v1".
                It is possible that "v_start" is not processed yet but "v1" has been processed:
                    - when the vertex goes to some end point and checks back.
                Therefore, an extra judgement is needed.
            --------------------------------------------------------------------------------------- */
            self.stack_push(v1 as usize);
        }

        // ::::::::::::::::::::::::::::::::::::::::
        //  Process this "v_start"
        // ::::::::::::::::::::::::::::::::::::::::
        the_p = v_start / 4; 
        op = self.op_string[the_p];
        remainder = (op % 4) as usize;

        // -----------------------------------------------------------------------------------------
        //  If this is a (valid) bond opeator, we flip it and the corresponding qudits
        // -----------------------------------------------------------------------------------------
        if remainder >= 2 {
            if self.flip == FLIPPED {
                self.op_string[the_p] = flip_operator!(self.op_string[the_p]);  
            }
            
        }

        // --------------------------------------------------------------------------
        //  Link to other vertices if "v_start" belongs to an off-diag site operator
        // --------------------------------------------------------------------------
        else {
            let v2: usize = v_start ^ 2;    // the neighbor
            let v3: usize = v_start ^ 1;    // the back
            let v4: usize = v3 ^ 2;         // the neighbor of the back

            if self.link_to_valid_dual_cluster_leg(v2) >= 0 {
                self.stack_push(v2);
            }

            if self.link_to_valid_dual_cluster_leg(v3) >= 0 {
                self.stack_push(v3);
            }

            if self.link_to_valid_dual_cluster_leg(v4) >= 0 {
                self.stack_push(v4);
            }
        }

        // Mark "v_start" as processed
        self.vertex_list[v_start] = self.flip;
    }
}

