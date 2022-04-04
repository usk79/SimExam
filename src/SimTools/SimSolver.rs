use std::collections::HashMap;

extern crate nalgebra as na;
use na::{U2, U3, Dynamic, ArrayStorage, VecStorage, Matrix, OMatrix, DMatrix};

pub mod simmodel;
use simmodel::{*};


#[derive(Debug)]
pub struct Simulator<T> 
where T: Model
{
    simtime: f64, // シミュレーション時間
    delta_t: f64, // 刻み幅Δt
    solvertype: SolverType, // 計算手法
    model: T,
    simstorage: HashMap<String, Vec<f64>>, // 計測する信号名とデータ配列の組み合わせ
}

impl<T> Simulator<T> 
where T: Model
{
    pub fn new(simtime: f64, delta_t: f64, solvertype: SolverType, model: T) -> Self {
        let stateinfo = model.get_state_info();
        let mut storage = HashMap::<String, Vec<f64>>::new();
        let storage_size = (simtime / delta_t + 0.5) as usize + 1;

        // シミュレーション結果を保存する領域を確保
        for e in stateinfo.iter() {
            storage.insert(e.to_string(), Vec::<f64>::with_capacity(storage_size));
        }

        Self {
            simtime: simtime,
            delta_t: delta_t,
            solvertype: solvertype,
            model: model,
            simstorage: storage,
        }
    }

    pub fn run_sim(&mut self) {
        let simsize = (self.simtime / self.delta_t + 0.5) as usize + 1;

        for idx in 0..simsize {
            
            self.model.calc_nextstate(self.delta_t, &self.solvertype);
            let vector = self.model.get_state().column(0);
            println!("{}, {}", vector[0], vector[1]);
        }
    }
}