use crate::tfim::TFIModel;

impl TFIModel {
    #[inline]
    pub fn rand_prob(&mut self) -> f64 {
        (self.rng.next() % u32::MAX) as f64 / (u32::MAX as f64)
    }

    #[inline]
    pub fn rand_bond(&mut self) -> usize {
        self.rng.next() as usize % self.num_bonds
    }

    #[inline]
    pub fn rand_site(&mut self) -> usize { 
        self.rng.next() as usize % self.num_sites 
    }

    #[inline]
    pub fn rand_qudit(&mut self) -> u8 { 
        (self.rng.next() as usize % 4) as u8 
    }
}