use std::io::{self, Write};
use std::fs::File;

use crate::types::{BasisState, LinearProgram, SimplexResult, SolveStatus};

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum LogLevel {
    Silent,
    Info,
    Debug,
}

enum Stream {
    Null,
    Stdout,
    File(File),
}

pub struct SimplexLogger {
    level: LogLevel,
    stream: Stream,
}

impl SimplexLogger {
    pub fn silent() -> Self {
        SimplexLogger { level: LogLevel::Silent, stream: Stream::Null }
    }

    pub fn new(level: LogLevel) -> Self {
        SimplexLogger { level, stream: Stream::Stdout }
    }

    pub fn to_file(path: &str, level: LogLevel) -> io::Result<Self> {
        let f = File::create(path)?;
        Ok(SimplexLogger { level, stream: Stream::File(f) })
    }

    fn write_line(&mut self, s: &str) {
        match &mut self.stream {
            Stream::Null => {}
            Stream::Stdout => println!("{}", s),
            Stream::File(f) => { let _ = writeln!(f, "{}", s); }
        }
    }

    pub fn log_problem(&mut self, lp: &LinearProgram) {
        if self.level == LogLevel::Silent { return; }
        self.write_line("=======================");
        self.write_line(" Revised Simplex Solver");
        self.write_line("=======================");
        self.write_line(&format!("A =\n{}", lp.a));
        self.write_line(&format!("b = {}", lp.b.transpose()));
        self.write_line(&format!("c = {}", lp.c.transpose()));
        self.write_line("");
    }

    pub fn log_phase(&mut self, phase: u8) {
        if self.level == LogLevel::Silent { return; }
        let header = if phase == 1 {
            "Phase 1 — Finding initial BFS"
        } else {
            "Phase 2 — Optimizing"
        };
        self.write_line(header);
        self.write_line(&"─".repeat(header.chars().count()));
    }

    pub fn log_iteration(&mut self, it: usize, basis: &BasisState, xb: &[f64], z: f64) {
        if self.level != LogLevel::Debug { return; }
        self.write_line(&format!("  iter {}:", it));
        self.write_line(&format!("    xb       = {:?}", xb));
        self.write_line(&format!("    basic    = {:?}", basis.basic));
        self.write_line(&format!("    nonbasic = {:?}", basis.nonbasic));
        self.write_line(&format!("    z        = {}", z));
        self.write_line("");
    }

    pub fn log_result(&mut self, result: &SimplexResult) {
        if self.level == LogLevel::Silent { return; }
        self.write_line("");
        let status_str = match result.status {
            SolveStatus::Optimal => "OPTIMAL",
            SolveStatus::Unbounded => "UNBOUNDED",
            SolveStatus::Infeasible => "INFEASIBLE",
        };
        self.write_line(&format!("Result: {}", status_str));
        self.write_line(&format!("  x          = {}", result.x.transpose()));
        self.write_line(&format!("  z          = {}", result.z));
        self.write_line(&format!("  iterations = {}", result.iterations));
    }
}
