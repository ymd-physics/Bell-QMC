use crate::tfim::TFIModel;
pub mod diagonal_update;
pub mod make_vertex_list;
pub mod cluster_update;
pub mod bond_cluster_update;
use crate::aux::{EMPTY, NULL_QUDIT, NULL_OP};

impl TFIModel {
    pub fn mc_thermalizing(&mut self) {
        self.diag_update();  
        self.cluster_update(); 
        self.refresh_left_right_qudits();    
        self.bond_cluster_update();
        self.adjust_m(); 
    }

    pub fn mc_sampling(&mut self) {
        self.diag_update_with_measure();   
        self.cluster_update(); 
        self.refresh_left_right_qudits();  
        self.bond_cluster_update();
    }

    fn adjust_m(&mut self) {
        let new_m = self.n + self.n / 3;

        if self.m < new_m {
            self.op_string.extend(vec![NULL_OP; new_m - self.m]);
            self.m = new_m;
            self.vertex_list = vec![EMPTY; 4 * new_m];
            self.stack = vec![0; 8 * new_m];
            self.left_qudits = vec![NULL_QUDIT; new_m];
            self.right_qudits = vec![NULL_QUDIT; new_m];     // "7" or "0b111" for invaild qudit
        }
    }
}