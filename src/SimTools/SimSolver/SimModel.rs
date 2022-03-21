extern crate nalgebra as na;
use na::{U2, U3, Dynamic, ArrayStorage, VecStorage, Matrix, OMatrix, DMatrix};

pub trait Model {
    fn slopefunc(&self) -> DMatrix<f64>;
}

#[derive(Debug, Clone)]
pub struct SpaceStateModel {
    mat_a: DMatrix<f64>,    // 状態遷移行列A
    mat_b: DMatrix<f64>,    // 入力行列B
    mat_c: DMatrix<f64>,    // 観測行列C
    state_dim: usize,       // 状態次数
    input_dim: usize,       // 入力次数
    output_dim: usize,      // 出力次数
    x: DMatrix<f64>,        // 状態ベクトル
    u: DMatrix<f64>,        // 入力ベクトル
}

impl SpaceStateModel {
    pub fn new(sdim: usize, idim: usize, odim: usize) -> Self {
        SpaceStateModel {
            mat_a: DMatrix::from_element(sdim, sdim, 0.0),
            mat_b: DMatrix::from_element(sdim, idim, 0.0),
            mat_c: DMatrix::from_element(odim, sdim, 0.0),
            x: DMatrix::from_element(sdim, 1, 0.0),
            u: DMatrix::from_element(idim, 1, 0.0),
            state_dim: sdim,
            input_dim: idim,
            output_dim: odim,
        }
    }

    pub fn init_state(&mut self, init_state: &[f64]) -> Result<(), &str> {

        if init_state.len() != self.state_dim {
            return Err("状態変数のサイズが違います。");
        }

        for (i, elem) in init_state.iter().enumerate() {
            self.x[i] = *elem;
        }

        Ok(())
    }

    pub fn set_mat_a(&mut self, mat_a: &[f64]) -> Result<(), &str> {

        if mat_a.len() != self.state_dim * self.state_dim {
            return Err("A行列のサイズが違います。");
        }

        for (i, elem) in mat_a.iter().enumerate() {
            self.mat_a[(i / self.state_dim, i % self.state_dim)] = *elem;
        }

        Ok(())
    }

    pub fn set_mat_b(&mut self, mat_b: &[f64]) -> Result<(), &str> {
        
        if mat_b.len() != self.state_dim * self.input_dim {
            return Err("B行列のサイズが違います。");
        }

        for (i, elem) in mat_b.iter().enumerate() {
            self.mat_b[(i / self.input_dim, i % self.input_dim)] = *elem;
        }

        Ok(())
    }

    pub fn set_mat_c(&mut self, mat_c: &[f64]) -> Result<(), &str> {

        if mat_c.len() != self.output_dim * self.state_dim {
            return Err("C行列のサイズが違います。");
        }

        for (i, elem) in mat_c.iter().enumerate() {
            self.mat_c[(i / self.state_dim, i % self.state_dim)] = *elem;
        }

        Ok(())
    }

}

impl Model for SpaceStateModel {
    fn slopefunc(&self) -> DMatrix<f64> {
        &self.mat_a * &self.x + &self.mat_b * &self.u
    }
}