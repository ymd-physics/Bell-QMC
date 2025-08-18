use std::time::Instant;
use std::fs::{File, OpenOptions};

pub const NULL_OP: i32 = -1;
pub const NULL_QUDIT: u8 = 7;

pub const FLIPPED: i32 = -2;
pub const NOT_FLIPPED: i32 = -1;

pub const EMPTY: i32 = -1;

pub const FREE_SPIN: i32 = -1;

#[macro_export]
macro_rules! get_rx {
    ($qudit: expr) => {
        $qudit & 1
    };
}

#[macro_export]
macro_rules! get_rz {
    ($qudit: expr) => {
        ($qudit >> 1) & 1
    };
}

#[macro_export]
macro_rules! go_through {
    ($leg: expr) => {
        $leg ^ 0b10
    };
}

#[macro_export]
macro_rules! move_to_neighbor {
    ($leg: expr) => {
        $leg ^ 0b01
    };
}

#[macro_export]
macro_rules! flip_rz {
    ($qudit: expr) => {
        $qudit ^ 0b10
    };
}

#[macro_export]
macro_rules! flip_rx {
    ($qudit: expr) => {
        $qudit ^ 0b01
    };
}


#[macro_export]
macro_rules! flip_operator {
    ($op: expr) => {
        $op ^ 1
    };
}

pub fn report_time(start_time: Instant) {
    let elapsed = start_time.elapsed().as_secs() as u32;
    let hours = elapsed / 3600;
    let minutes = (elapsed % 3600) / 60;
    let seconds = elapsed % 60;
    println!("■ Runtime: {}h-{}min-{}s", hours, minutes, seconds);
}

pub fn print_horizontal_line(len: usize, marker: &str) {
    for _ in 0..len { print!("{}", marker); } print!("\n");
}

#[macro_export]
macro_rules! dual_go_through {
    ($leg: expr) => {
        $leg ^ 0b01
    };
}

pub fn create_new_file(file_path: String) -> File {
    OpenOptions::new().write(true).append(true).create(true).open(file_path).unwrap()
}