use crate::tfim::TFIModel;
/***********************************************************
 *      For s = (s^z, s^x), the Pauli matrix is 
 *              00 ~ I
 *              01 ~ X
 *              10 ~ Z
 *              11 ~ Y
 **********************************************************/
impl TFIModel {
    // ==================================================
    //  Pauli operators
    // ==================================================
    const BASE: f64 = -1.0;

    // #[inline]
    // fn measure_sigma(&self, sigma: u8, s: usize) -> f64 {
    //     assert!(sigma < 4, "Invalid Pauli matrix specified.");
    //     let s_z: u8 = (sigma >> 1) & 1;
    //     let s_x: u8 = sigma & 1;
    //     let r_z: u8 = (self.qudits[s] >> 1 ) & 1;
    //     let r_x: u8 = self.qudits[s] & 1;
    //     Self::BASE.powf((s_x * r_z - s_z * r_x) as f64)
    // }

    #[inline]
    fn measure_x(&self, s: usize) -> f64 {
        let r_z: u8 = (self.qudits[s] >> 1 ) & 1;
        Self::BASE.powf(r_z as f64)
    }

    // #[inline]
    // fn measure_y(&self, s: usize) -> f64 {
    //     let r_z: u8 = (self.qudits[s] >> 1 ) & 1;
    //     let r_x: u8 = self.qudits[s] & 1;
    //     Self::BASE.powf((r_z - r_x) as f64)
    // }

    #[inline]
    fn measure_z(&self, s: usize) -> f64 {
        let r_x: u8 = self.qudits[s] & 1;
        Self::BASE.powf(-1.0 * r_x as f64)
    }

    // ==================================================
    //  Purity
    // ==================================================
    #[inline]
    fn measure_swap(&self, s: usize) -> f64 {
        if self.qudits[s] == 0b11 { -1.0 } else { 1.0 }
    }

    fn measure_purity(&self, swapped_region: &[usize]) -> f64 {
        swapped_region.iter().map(|&s| self.measure_swap(s)).product()
    }

    // ==================================================
    //  Pauli correlations
    // ==================================================
    #[inline]
    fn get_xx_corr2_(&self, s_i: usize, s_j: usize) -> f64 {
        self.measure_x(s_i) * self.measure_x(s_j)
    }

    #[inline]
    fn get_zz_corr_2(&self, s_i: usize, s_j: usize) -> f64 {
        self.measure_z(s_i) * self.measure_z(s_j)
    }    

    
}

impl TFIModel {
    pub fn ini_measure(&mut self) {
        self.purity = 0.0;
        self.partial_purity = 0.0;

        for s in 0..self.num_sites {
            self.zz_corr_2[s] = 0.0;
            self.xx_corr_2[s] = 0.0;
        }
    }

    pub fn measure(&mut self) {
        self.purity += self.measure_purity(&self.system);
        self.partial_purity += self.measure_purity(&self.subsystem);

        for s in 0..self.num_sites {
            self.zz_corr_2[s] += self.get_zz_corr_2(0, s);
            self.xx_corr_2[s] += self.get_xx_corr2_(0, s);
        }
    }

    pub fn statisticize(&mut self, num_samples: f64) {
        self.purity /= num_samples;
        self.partial_purity /= num_samples;

        for s in 0..self.num_sites {
            self.zz_corr_2[s] /= num_samples;
            self.xx_corr_2[s] /= num_samples; 
        }
    }
}