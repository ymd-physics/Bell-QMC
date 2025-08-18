use crate::tfim::TFIModel;

impl TFIModel {
    #[inline]
    pub fn stack_initialize(&mut self) {
        self.top = -1;
    }

    #[inline]
    pub fn stack_push(&mut self, x: usize) {
        self.top += 1;
        self.stack[self.top as usize] = x;
    }

    #[inline]
    pub fn stack_pop(&mut self) -> usize {
        let top_value = self.stack[self.top as usize];
        self.top -= 1;
        top_value
    }
}