use crate::tfim::TFIModel;
use crate::aux::{FLIPPED, NOT_FLIPPED, EMPTY, FREE_SPIN};
use crate::{flip_operator, flip_rx, get_rz, go_through};

impl TFIModel {
    fn link_to_valid_cluster_leg(&mut self, v: usize) -> i32 {
        // --------------------------------------------------
        //  It stops at 
        //      (i) a flippable site operator
        //      (ii) an off-diagonal bond operator
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
                
                // an off-diagonal bond operator is always okay
                // A diag-site operator that dose not satisfy the constraint is not valid
                if (remainder1 == 2) || ((remainder1 == 0) && (get_rz!(self.left_qudits[the_p1]) == 1)) {
                        v0 = go_through!(v1); 
                    }

                else {
                    break;
                }
            }

            v1 as i32 
        }
    }
    
    pub fn cluster_update(&mut self) {
        self.make_vertex_list();

        let mut op: i32;
        let mut the_p: usize;
        let mut remainder: usize;

        let mut v0: usize;
        self.stack_initialize();        // for growing the cluster

        for v in (0..4 * self.m).step_by(2) {
            if self.vertex_list[v] < 0 { 
                continue; 
            }

            // --------------------------------------------------------------------
            // We must ensure it is a valid leg, i.e. it belongs to 
            //  (i) a flippable site operator (ii) an off-diagonal bond operator
            // --------------------------------------------------------------------
            the_p = v / 4;
            op = self.op_string[the_p];
            remainder = (op % 4) as usize;

            // --------------------------------------------------------------------
            // an off-diagonal site operator is always valid
            // But for a diag-site operator, we skip it if r^z = 1
            // --------------------------------------------------------------------
            if remainder < 2 {
                if (remainder == 1) || (get_rz!(self.left_qudits[the_p]) == 0) {
                    self.flip = if self.rand_prob() > 0.5 { FLIPPED } else { NOT_FLIPPED };
                    self.stack_push(v);

                    // --------------------------------------------
                    //  Grow the cluster starting with "v"
                    // --------------------------------------------
                    loop {
                        if self.top == EMPTY { break; }
                        self.make_cluster();
                    }
                }
            }            
        }

        // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        //  Update the qudits at zero time after the updates of site operators
        // :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        for s in 0..self.num_sites {
            if self.v_first[s] != FREE_SPIN {
                let v = self.v_first[s] as usize;  
                v0 = v;     // this is the start point

                loop {
                    if self.vertex_list[v0] < 0 {
                        if self.vertex_list[v0] == FLIPPED {
                            self.qudits[s] = flip_rx!(self.qudits[s]);
                        }
                        break;
                    }
                    
                    else {
                        v0 = self.vertex_list[go_through!(v0)] as usize;    
                        if v0 == v {
                            break;
                        }
                    }
                }
            }

            else {
                // Randomly change those isolated qudits (without changing the sector)
                if self.rand_prob() > 0.5 { self.qudits[s] = flip_rx!(self.qudits[s]); }
            }
        }
    }

    fn make_cluster(&mut self) {
        let the_p: usize;
        let op: i32;
        let remainder: usize;

        let v_start: usize = self.stack_pop();
        let v1: i32 = self.link_to_valid_cluster_leg(v_start);

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

        // -----------------------------------------------------
        //  If this is a (valid) site opeator, and we flip it 
        // -----------------------------------------------------
        if remainder < 2 {
            if self.flip == FLIPPED {
                self.op_string[the_p] = flip_operator!(self.op_string[the_p]);
            }     
        }
        
        // --------------------------------------------------------------------------
        //  Link to other vertices if "v_start" belongs to an off-diag bond operator
        // --------------------------------------------------------------------------
        else {
            let v2: usize = v_start ^ 1;    // the neighbor
            let v3: usize = v_start ^ 2;    // the back
            let v4: usize = v3 ^ 1;         // the neighbor of the back

            if self.link_to_valid_cluster_leg(v2) >= 0 {
                self.stack_push(v2);
            }

            if self.link_to_valid_cluster_leg(v3) >= 0 {
                self.stack_push(v3);
            }

            if self.link_to_valid_cluster_leg(v4) >= 0 {
                self.stack_push(v4);
            }
        }

        // ---------------------------------
        //  Mark "v_start" as processed
        // ---------------------------------
        self.vertex_list[v_start] = self.flip;
    }

}