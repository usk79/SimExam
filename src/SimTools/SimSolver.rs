use std::fs;
use std::fs::File;
use std::io::{self, Write, BufWriter};

use std::collections::HashMap;

extern crate nalgebra as na;
use na::{U2, U3, Dynamic, ArrayStorage, VecStorage, Matrix, OMatrix, DMatrix};

use plotters::prelude::*;


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

    pub fn timeplot(&self, dirname : &str, pltsize: (u32, u32)) {
        match fs::create_dir(dirname) { 
            Err(e) => println!("! {:?}", e.kind()),
            Ok(_) => {},
        }

        let signalinfo = self.model.get_signals_info();
        for signal in signalinfo.iter() {
            self.timeplot_subfn(dirname, &signal, pltsize);
        }

    }

    fn timeplot_subfn(&self, dirname: &str, signal: &String, pltsize: (u32, u32)) {
        let filename = format!("./{}/{}.svg", dirname, signal);

        let plt = BitMapBackend::new(&filename, pltsize).into_drawing_area();
        plt.fill(&WHITE).unwrap();
    
        let font = ("sans-serif", 20);

        let timeaxis: &Vec<f64> = self.simstorage.get("time").unwrap();
        let valuaxis: &Vec<f64> = self.simstorage.get(signal).unwrap();
        let (y_min, y_max) = valuaxis.iter()
                         .fold(
                           (0.0/0.0, 0.0/0.0),
                           |(m,n), v| (v.min(m), v.max(n))
                          ); // f64はNaNがあるためordが実装されていない。min, maxを使うための工夫が必要⇒https://qiita.com/lo48576/items/343ca40a03c3b86b67cb

        let xrange = 0.0..self.simtime; 
        let yrange = y_min..y_max;
      
        let mut chart = ChartBuilder::on(&plt)
          .caption(&signal, font.into_font()) // キャプションのフォントやサイズ
          .margin(10)                         // 上下左右全ての余白
          .x_label_area_size(16)              // x軸ラベル部分の余白
          .y_label_area_size(42)              // y軸ラベル部分の余白
          .build_cartesian_2d(                // x軸とy軸の数値の範囲を指定する
            xrange,                           // x軸の範囲
            yrange)                           // y軸の範囲
          .unwrap();
    
        // x軸y軸、グリッド線などを描画
        chart.configure_mesh().draw().unwrap();

        let line_series = LineSeries::new(
            timeaxis.iter()
                    .zip(valuaxis.iter())
                    .map(|(x, y)| (*x, *y)),
                &RED);

        chart.draw_series(line_series).unwrap();

    }
}