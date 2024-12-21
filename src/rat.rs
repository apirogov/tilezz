use crate::snake::Snake;
use crate::zzbase::ZZNum;

pub struct Rat<T: ZZNum> {
    pub snake: Snake<T>,
}

impl<T: ZZNum> Rat<T> {
    pub fn new(snake: Snake<T>) -> Self {
        assert!(snake.is_rat());
        Self { snake: snake }
    }
}
