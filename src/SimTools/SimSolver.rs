extern crate nalgebra as na;
use na::{U2, U3, Dynamic, ArrayStorage, VecStorage, Matrix, OMatrix, DMatrix};

pub mod SimModel;
use SimModel::{*};

#[derive(Debug)]
pub enum SolverType {
    Euler,
    RungeKutta,
}

#[derive(Debug)]
pub struct Simulator<T> 
where T: Model
{
    simtime: f64, // シミュレーション時間
    delta_t: f64, // 刻み幅Δt
    solvertype: SolverType, // 計算手法
    models: Vec<T>
}

impl<T> Simulator<T> 
where T: Model
{
    pub fn new(simtime: f64, delta_t: f64, solvertype: SolverType) -> Self {
        Self {
            simtime: simtime,
            delta_t: delta_t,
            solvertype: solvertype,
            models: Vec::new(),
        }
    }

    pub fn run_sim(&self) {
        match &self.solvertype {
            SolverType::Euler       => self.euler_method(),
            SolverType::RungeKutta  => self.runge_kutta_method(),
        }
    }

    fn euler_method(&self) {
        println!("euler");
    }

    fn runge_kutta_method(&self) {
        println!("runge_kutta");
    }
}