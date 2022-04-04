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
        storage.insert("time".to_string(), Vec::<f64>::with_capacity(storage_size)); // 時刻用
        storage.get_mut("time").unwrap().push(0.0); // simstorageに時刻0のデータを格納

        let state = model.get_state();
        for (i, e) in stateinfo.iter().enumerate() {
            storage.insert(e.to_string(), Vec::<f64>::with_capacity(storage_size));
            storage.get_mut(&e.to_string()).unwrap().push(state[i]);
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
        let stateinfo = self.model.get_state_info();

        for idx in 1..simsize {
            self.simstorage.get_mut("time").unwrap().push(idx as f64 * self.delta_t); // 時刻の記録

            self.model.calc_nextstate(self.delta_t, &self.solvertype); // 1ステップ進める
            let state = self.model.get_state(); // 現在の状態を取得する 
            for (i, e) in stateinfo.iter().enumerate() {
                self.simstorage.get_mut(&e.to_string()).unwrap().push(state[i]); // 計算結果を記録
            }
        }

        println!("{:?}", self.simstorage);
    }
}