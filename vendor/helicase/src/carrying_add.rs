#[cfg(target_arch = "x86_64")]
mod x86 {
    #[derive(Default, Clone, Copy)]
    pub struct Carry(u8);

    impl Carry {
        #[inline(always)]
        pub const fn new(carry: bool) -> Self {
            Self(carry as u8)
        }

        #[allow(unused_unsafe)]
        #[inline(always)]
        pub fn add(&mut self, lhs: u64, rhs: u64) -> u64 {
            use core::arch::x86_64::_addcarry_u64;
            let mut res = 0;
            self.0 = unsafe { _addcarry_u64(self.0, lhs, rhs, &mut res) };
            res
        }
    }
}

#[cfg(target_arch = "x86_64")]
pub use x86::*;

#[cfg(not(target_arch = "x86_64"))]
mod fallback {
    #[derive(Default, Clone, Copy)]
    pub struct Carry(bool);

    impl Carry {
        #[inline(always)]
        pub const fn new(carry: bool) -> Self {
            Self(carry)
        }

        // Adapted from https://github.com/rust-lang/rust/blob/master/library/core/src/num/uint_macros.rs#L2428
        #[inline(always)]
        pub const fn add(&mut self, lhs: u64, rhs: u64) -> u64 {
            let (a, c1) = lhs.overflowing_add(rhs);
            let (b, c2) = a.overflowing_add(self.0 as u64);
            self.0 = c1 | c2;
            b
        }
    }
}

#[cfg(not(target_arch = "x86_64"))]
pub use fallback::*;
