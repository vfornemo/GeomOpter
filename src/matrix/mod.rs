use std::{fmt::Debug, ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign}, path::Iter};


/// Column major matrix
/// 
/// Example:
/// size: [2,3]
/// data: [1,2,3,4,5,6]
/// 
/// 1 3 5
/// 
/// 2 4 6
/// 
#[derive(Debug, Clone)]
pub struct MatFull<T> {
    pub size: [usize;2],
    pub data: Vec<T>,
}


impl <T> MatFull<T> 
    where T: Copy + Clone + Debug + PartialEq
{

    pub fn new(size: [usize;2]) -> MatFull<T> {
        MatFull {
            size: [0,0],
            data: vec![],
        }
    }

    pub fn from_vec(size: [usize;2], data: Vec<T>) -> MatFull<T> {
        MatFull {
            size,
            data,
        }
    }

    pub fn reshape(&mut self, size: [usize;2]) {
        self.size = size;
    }

    pub fn get(&self, i: usize, j: usize) -> &T {
        &self.data[j*self.size[0] + i]
    }

    pub fn get_mut(&mut self, i: usize, j: usize) -> &mut T {
        &mut self.data[j*self.size[0] + i]
    }

    pub fn iter_row(&self, i: usize) -> impl Iterator<Item = T> + '_  {
        self.data.iter().enumerate().filter(move |(idx,_x)| idx % self.size[0] == i)
                .map(|(_idx,x)| *x)
    }

    pub fn iter_col(&self, j: usize) -> std::slice::Iter<T> {
        self.data[j*self.size[0]..(j+1)*self.size[0]].iter()
    }

    pub fn transpose(&self) -> MatFull<T> {
        let mut data = Vec::new();
        for i in 0..self.size[0] {
            for j in 0..self.size[1] {
                data.push(self.get(i, j).clone());
            }
        }
        MatFull {
            size: [self.size[1], self.size[0]],
            data,
        }
    }


    pub fn formated_output(&self) {
        let mut output = String::new();
        for i in 0..self.size[0] {
            for j in 0..self.size[1] {
                output.push_str(&format!("{:?}", self.get(i, j)));
                output.push_str(" ");
            }
            output.push_str("\n");
        }
        
        println!("{:?}", output);
    }
    

    pub fn mat_mul(&self, other: &MatFull<T>) -> MatFull<T> 
        where T: Add<Output = T> + Mul<Output = T> + Copy + Clone + Debug + PartialEq + std::iter::Sum
    {
        if self.size[1] != other.size[0] {
            panic!("Matrix size mismatch: {:?} != {:?}", self.size, other.size);
        }
        let size_new = [self.size[0], other.size[1]];
        let mut data = Vec::new();
        for j in 0..other.size[1] {
            for i in 0..self.size[0] {
                data.push(self.iter_row(i).zip(other.iter_col(j)).map(|(x,y)| x* *y).sum())
            }
        }

        MatFull {
            size: size_new,
            data,
        }
    }



    
}

impl <T> AddAssign for MatFull<T> 
    where T: AddAssign + Add + Copy + Clone + Debug + PartialEq
{
    fn add_assign(&mut self, other: Self) {
        if self.size != other.size {
            panic!("Matrix size mismatch: {:?} != {:?}", self.size, other.size);
        }
        for i in 0..self.data.len() {
            self.data[i] += other.data[i];
        }
    }
}


impl <T> SubAssign for MatFull<T> 
    where T: SubAssign + Sub + Copy + Clone + Debug + PartialEq
{
    fn sub_assign(&mut self, other: Self) {
        if self.size != other.size {
            panic!("Matrix size mismatch: {:?} != {:?}", self.size, other.size);
        }
        for i in 0..self.data.len() {
            self.data[i] -= other.data[i];
        }
    }
}

impl MulAssign<f64> for MatFull<f64> 
{
    fn mul_assign(&mut self, a: f64) {
        for i in 0..self.data.len() {
            self.data[i] *= a;
        }

    }
}






#[test]
fn test_matrix() {
    let size = [2,3];
    let data = vec![1,2,3,4,5,6];
    let mut mat = MatFull::from_vec(size, data);
    assert_eq!(*mat.get(0,0), 1);
    assert_eq!(*mat.get(1,2), 6);
    *mat.get_mut(0,0) = 10;
    assert_eq!(*mat.get(0,0), 10);
    mat.reshape([3,2]);
    assert_eq!(*mat.get(0,0), 10);
    assert_eq!(*mat.get(2,1), 6);
    let mat_t = mat.transpose();
    assert_eq!(mat_t.data, vec![10,4,2,5,3,6]);
    mat_t.formated_output();    
}

#[test]
fn test_mat_ops() {
    let size = [3,2];
    let data = vec![1,2,3,4,5,6];
    let mut mat1 = MatFull::from_vec(size, data);
    let size = [2,2];
    let data = vec![1,2, 3, 4];
    let mat2 = MatFull::from_vec(size, data);
    let mat3 = mat1.mat_mul(&mat2);
    mat3.formated_output();
    let mat3 = MatFull::from_vec([3,2], vec![1,2,3,4,5,6]);
    mat1 += mat3.clone();
    mat1.formated_output();
    mat1 -= mat3.clone();
    mat1.formated_output();
    let mut mat4 = MatFull::from_vec([3,2], vec![1.0,2.0,3.0,4.0,5.0,6.0]);
    mat4 *= 2.0;
    mat4.formated_output();
}