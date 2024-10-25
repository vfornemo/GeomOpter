use std::fmt::Debug;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Sub, SubAssign};


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

    pub fn new() -> MatFull<T> {
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

    pub fn get_ij(&self, i: usize, j: usize) -> &T {
        &self.data[j*self.size[0] + i]
    }

    pub fn get_mut_ij(&mut self, i: usize, j: usize) -> &mut T {
        &mut self.data[j*self.size[0] + i]
    }

    pub fn iter_row(&self, i: usize) -> impl Iterator<Item = T> + '_  {
        self.data.iter().enumerate().filter(move |(idx,_x)| idx % self.size[0] == i)
                .map(|(_idx,x)| *x)
    }

    pub fn iter_col(&self, j: usize) -> std::slice::Iter<T> {
        self.data[j*self.size[0]..(j+1)*self.size[0]].iter()
    }

    pub fn iter_col_mut(&mut self, j: usize) -> std::slice::IterMut<T> {
        self.data[j*self.size[0]..(j+1)*self.size[0]].iter_mut()
    }

    pub fn iter_col_range(&self, j: usize, start: usize, end: usize) -> std::slice::Iter<T> {
        self.data[j*self.size[0] + start..j*self.size[0] + end].iter()
    }

    pub fn iter_col_range_mut(&mut self, j: usize, start: usize, end: usize) -> std::slice::IterMut<T> {
        self.data[j*self.size[0] + start..j*self.size[0] + end].iter_mut()
    }

    pub fn transpose(&self) -> MatFull<T> {
        let mut data = Vec::new();
        for i in 0..self.size[0] {
            for j in 0..self.size[1] {
                data.push(self.get_ij(i, j).clone());
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
                output.push_str(&format!("{:?}", self.get_ij(i, j)));
                output.push_str(" ");
            }
            output.push_str("\n");
        }
        
        println!("{}", output);
    }
    
    // /// change to n,
    // pub fn mat_mul(&self, other: &MatFull<T>, op_a: char, op_b: char) -> MatFull<T> 
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

impl MatFull<f64> {

    /// Formated print for matrix
    /// Input:
    /// op: 't' for transpose, 'n' for normal
    /// precision: number of decimal places
    /// 
    pub fn formated_output_f64(&self, op: char, precision: usize) {
        let mut output = String::new();
        if op == 't' {
            for i in 0..self.size[1] {
                for j in 0..self.size[0] {
                    output.push_str(&format!("{:-11.*}", precision, self.get_ij(j, i)));
                    output.push_str(" ");
                }
                output.push_str("\n");
            }

        } else {
            for i in 0..self.size[0] {
                for j in 0..self.size[1] {
                    output.push_str(&format!("{:-11.*}", precision, self.get_ij(i, j)));
                    output.push_str(" ");
                }
                output.push_str("\n");
            }
        }



        
        println!("{}", output);

    }

}

impl <'a, 'b, T> Add for &'a MatFull<T> 
    where T: Add<Output = T> + Copy + Clone + Debug + PartialEq
{
    type Output = MatFull<T>;

    fn add(self, other: &MatFull<T>) -> MatFull<T> {
        if self.size != other.size {
            panic!("Matrix size mismatch: {:?} != {:?}", self.size, other.size);
        }
        let data = self.data.iter().zip(other.data.iter()).map(|(x,y)| (*x + *y)).collect();

        MatFull {
            size: self.size,
            data,
        }
    }
    
}

impl <T> Add for MatFull<T> 
    where T: Add<Output = T> + Copy + Clone + Debug + PartialEq
{
    type Output = MatFull<T>;
    /// T + U
    fn add(self, other: Self) -> Self::Output {
        if self.size != other.size {
            panic!("Matrix size mismatch: {:?} != {:?}", self.size, other.size);
        }
        let data = self.data.iter().zip(other.data.iter()).map(|(x,y)| (*x + *y)).collect();

        MatFull {
            size: self.size,
            data,
        }
    }
    
}

impl <T> AddAssign<&MatFull<T>> for MatFull<T> 
    where T: AddAssign + Add + Copy + Clone + Debug + PartialEq
{
    /// &T + &U
    fn add_assign(&mut self, other: &MatFull<T>) {
        if self.size != other.size {
            panic!("Matrix size mismatch: {:?} != {:?}", self.size, other.size);
        }
        self.data.iter_mut().zip(other.data.iter()).for_each(|(x,y)| *x += *y);
    }
}

// impl <T> AddAssign for MatFull<T> 
//     where T: AddAssign + Add + Copy + Clone + Debug + PartialEq
// {
//     /// &T + U
//     fn add_assign(&mut self, other: MatFull<T>) {
//         if self.size != other.size {
//             panic!("Matrix size mismatch: {:?} != {:?}", self.size, other.size);
//         }
//         self.data.iter_mut().zip(other.data.iter()).for_each(|(x,y)| *x += *y);
//     }
// }

//


impl <T> Sub for MatFull<T> 
    where T: Sub<Output = T> + Copy + Clone + Debug + PartialEq
{
    type Output = MatFull<T>;
    /// T - U
    fn sub(self, other: Self) -> Self::Output {
        if self.size != other.size {
            panic!("Matrix size mismatch: {:?} != {:?}", self.size, other.size);
        }
        let data = self.data.iter().zip(other.data.iter()).map(|(x,y)| (*x - *y)).collect();
        MatFull {
            size: self.size,
            data,
        }
    }
    
}

// impl <T> SubAssign for MatFull<T> 
//     where T: SubAssign + Sub + Copy + Clone + Debug + PartialEq
// {
//     /// &T - U
//     fn sub_assign(&mut self, other: Self) {
//         if self.size != other.size {
//             panic!("Matrix size mismatch: {:?} != {:?}", self.size, other.size);
//         }
//         self.data.iter_mut().zip(other.data.iter()).for_each(|(x,y)| *x -= *y);
//     }
// }

impl <T> SubAssign<&MatFull<T>> for MatFull<T> 
    where T: SubAssign + Sub + Copy + Clone + Debug + PartialEq
{
    /// &T - &U
    fn sub_assign(&mut self, other: &MatFull<T>) {
        if self.size != other.size {
            panic!("Matrix size mismatch: {:?} != {:?}", self.size, other.size);
        }
        self.data.iter_mut().zip(other.data.iter()).for_each(|(x,y)| *x -= *y);
    }
}

impl MulAssign<f64> for MatFull<f64> 
{
    fn mul_assign(&mut self, a: f64) {
        self.data.iter_mut().for_each(|x| *x *= a);

    }
}

impl Mul<f64> for MatFull<f64> 
{
    type Output = MatFull<f64>;

    fn mul(self, a: f64) -> Self::Output {
        let data = self.data.iter().map(|x| *x * a).collect();
        MatFull {
            size: self.size,
            data,
        }
    }
    
}

impl DivAssign<f64> for MatFull<f64> 
{
    fn div_assign(&mut self, a: f64) {
        self.data.iter_mut().for_each(|x| *x /= a);

    }
}

impl Div<f64> for MatFull<f64> 
{
    type Output = MatFull<f64>;

    fn div(self, a: f64) -> Self::Output {
        let data = self.data.iter().map(|x| *x / a).collect();
        MatFull {
            size: self.size,
            data,
        }
    }
    
}

#[test]
fn test_matrix() {
    let size = [2,3];
    let data = vec![1,2,3,4,5,6];
    let mut mat = MatFull::from_vec(size, data);
    assert_eq!(*mat.get_ij(0,0), 1);
    assert_eq!(*mat.get_ij(1,2), 6);
    *mat.get_mut_ij(0,0) = 10;
    assert_eq!(*mat.get_ij(0,0), 10);
    mat.reshape([3,2]);
    assert_eq!(*mat.get_ij(0,0), 10);
    assert_eq!(*mat.get_ij(2,1), 6);
    let mat_t = mat.transpose();
    assert_eq!(mat_t.data, vec![10,4,2,5,3,6]);
    mat_t.formated_output();    

    let datan = vec![1,2,3,4,5,6, 7, 8, 9];
    let matn = MatFull::from_vec([3,3], datan);

    let slice: Vec<i32> = matn.iter_col_range(1, 0, 2).map(|x| *x).collect();
    println!("slice {:?}", slice);
    
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
    mat1 += &mat3;
    mat1.formated_output();
    mat1 -= &mat3;
    mat1.formated_output();
    let mut mat4 = MatFull::from_vec([3,2], vec![1.0,2.0,3.0,4.0,5.0,6.0]);
    mat4 *= 2.0;
    mat4.formated_output();
}