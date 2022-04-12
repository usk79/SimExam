use std::fs::File;
use std::io::{self, Write, BufWriter};

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
    simsize: usize,
    solvertype: SolverType, // 計算手法
    model: T,
    simstorage: HashMap<String, Vec<f64>>, // 計測する信号名とデータ配列の組み合わせ
}

impl<T> Simulator<T> 
where T: Model
{
    pub fn new(simtime: f64, delta_t: f64, solvertype: SolverType, model: T) -> Self {
        let signalinfo = model.get_signals_info();
        let mut storage = HashMap::<String, Vec<f64>>::new();
        let storage_size = (simtime / delta_t + 0.5) as usize + 1;
        let simsize = (simtime / delta_t + 0.5) as usize + 1;

        // シミュレーション結果を保存する領域を確保
        storage.insert("time".to_string(), Vec::<f64>::with_capacity(storage_size)); // 時刻用
        storage.get_mut("time").unwrap().push(0.0); // simstorageに時刻0のデータを格納

        let signals = model.get_allsignals();
        for (i, e) in signalinfo.iter().enumerate() {
            storage.insert(e.to_string(), Vec::<f64>::with_capacity(storage_size));
            storage.get_mut(&e.to_string()).unwrap().push(signals[i]);
        }

        Self {
            simtime: simtime,
            delta_t: delta_t,
            simsize: simsize,
            solvertype: solvertype,
            model: model,
            simstorage: storage,
        }
    }

    pub fn run_sim(&mut self) {
        let signalinfo = self.model.get_signals_info();

        for idx in 1..self.simsize {
            self.simstorage.get_mut("time").unwrap().push(idx as f64 * self.delta_t); // 時刻の記録

            self.model.calc_nextstate(self.delta_t, &self.solvertype); // 1ステップ進める
            let signals = self.model.get_allsignals(); // 現在の状態を取得する 
            for (i, e) in signalinfo.iter().enumerate() {
                self.simstorage.get_mut(&e.to_string()).unwrap().push(signals[i]); // 計算結果を記録
            }
        }
    }

    pub fn export_sim(&self, filepath: &str) { // csv形式として吐き出す
        let mut file = BufWriter::new(File::create(filepath).unwrap());

        // 先頭の行（系列名を書き出す）
        //let keyary = self.simstorage.iter().map(|(k, _v)| k.to_string()).collect::<Vec<String>>(); // HashMapは順番が保持されない
        let mut seriesname = vec![String::from("time")];
        seriesname.append( &mut self.model.get_signals_info() );
        writeln!(file, "{}", seriesname.join(","));

        for idx in 0..self.simsize { 
            // 1行ずつ作成
            let line = seriesname.iter().map(|name| self.simstorage.get(name).unwrap()[idx].to_string())
                                        .collect::<Vec<String>>().join(",");
            
            writeln!(file, "{}", line).unwrap();
        }
        
    }
}