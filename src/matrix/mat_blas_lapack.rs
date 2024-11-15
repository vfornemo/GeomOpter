use blas::dgemm;
use lapack::dsyev;
use crate::matrix::MatFull;
extern crate rand;

/// BLAS dgemm 
fn _dgemm<'a>(mat_a: &'a MatFull<f64>, mat_b: &'a MatFull<f64>, mat_c: &'a mut MatFull<f64>, opa: char, opb: char, alpha: f64, beta: f64) {	

    let size_check = match (&opa, &opb) {
        ('N', 'N') => mat_a.size[1] == mat_b.size[0] && mat_a.size[0] == mat_c.size[0] && mat_b.size[1] == mat_c.size[1],
        ('n', 'n') => mat_a.size[1] == mat_b.size[0] && mat_a.size[0] == mat_c.size[0] && mat_b.size[1] == mat_c.size[1],
        ('N', 'T') => mat_a.size[1] == mat_b.size[1] && mat_a.size[0] == mat_c.size[0] && mat_b.size[0] == mat_c.size[1],
        ('n', 't') => mat_a.size[1] == mat_b.size[1] && mat_a.size[0] == mat_c.size[0] && mat_b.size[0] == mat_c.size[1],
        ('T', 'N') => mat_a.size[0] == mat_b.size[0] && mat_a.size[1] == mat_c.size[0] && mat_b.size[1] == mat_c.size[1],
        ('t', 'n') => mat_a.size[0] == mat_b.size[0] && mat_a.size[1] == mat_c.size[0] && mat_b.size[1] == mat_c.size[1],
        ('T', 'T') => mat_a.size[0] == mat_b.size[1] && mat_a.size[1] == mat_c.size[0] && mat_b.size[0] == mat_c.size[1],
        ('t', 't') => mat_a.size[0] == mat_b.size[1] && mat_a.size[1] == mat_c.size[0] && mat_b.size[0] == mat_c.size[1],
        _ => false,
    };

    if !size_check {
        panic!("Matrix size mismatch or invalid operation in dgemm. A: {:?}, B: {:?}, C: {:?}", mat_a.size, mat_b.size, mat_c.size);
    }

    let m = if opa == 'N' || opa == 'n' { mat_a.size[0] } else { mat_a.size[1] };
    let k = if opa == 'N' || opa == 'n' { mat_a.size[1] } else { mat_a.size[0] };
    let n = if opb == 'N' || opb == 'n' { mat_b.size[1] } else { mat_b.size[0] };
    let lda =  if opa == 'N' || opa == 'n' { m.max(1) } else { k.max(1) };
    let ldb =  if opb == 'N' || opb == 'n' { k.max(1) } else { n.max(1) };

    unsafe {
        dgemm(opa as u8,
              opb as u8, 
              m as i32,
              n as i32,
              k as i32,
              alpha,
              mat_a.data.as_ref(),
              lda as i32,
              mat_b.data.as_ref(),
              ldb as i32,
              beta,
              mat_c.data.as_mut(),
              mat_c.size[0] as i32);
    }
                

}

/// Matrix multiplication 
/// ```text
/// C := alpha*op( A )*op( B ) + beta*C, 
/// ```
/// opa: `'N'/'n'`(normal) or `'T'/'t'`(transpose) for A \
/// opb: `'N'/'n'`(normal) or `'T'/'t'`(transpose) for B \
pub fn mat_dgemm<'a>(mat_a: &'a MatFull<f64>, mat_b: &'a MatFull<f64>, opa: char, opb: char, alpha: f64, beta: f64) -> MatFull<f64> {
    let new_size = match (&opa, &opb) {
        ('N', 'N') => [mat_a.size[0], mat_b.size[1]],
        ('n', 'n') => [mat_a.size[0], mat_b.size[1]],
        ('N', 'T') => [mat_a.size[0], mat_b.size[0]],
        ('n', 't') => [mat_a.size[0], mat_b.size[0]],
        ('T', 'N') => [mat_a.size[1], mat_b.size[1]],
        ('t', 'n') => [mat_a.size[1], mat_b.size[1]],
        ('T', 'T') => [mat_a.size[1], mat_b.size[0]],
        ('t', 't') => [mat_a.size[1], mat_b.size[0]],
        _ => panic!("Invalid matrix operation"),
    };
    let mut mat_c = MatFull::new_with_value(new_size, 0.0);
    _dgemm(mat_a, mat_b, &mut mat_c, opa, opb, alpha, beta);

    mat_c
    
}


/// Calculate the eigenvalues and eigenvectors of a real symmetric matrix \
/// jobz: `'N'/'n'` (eigenvalues only) or `'V'/'v'` (eigenvalues and eigenvectors)
pub fn mat_dsyev<'a>(mat_a: &'a MatFull<f64>, jobz: char) -> (Option<MatFull<f64>>, Vec<f64>, i32) {
    if mat_a.size[0] != mat_a.size[1] {
        panic!("Error in _dsyev: not a symmetric matrix: {:?}", mat_a.size);
    }

    let mut a = mat_a.data.iter().map(|x| *x).collect::<Vec<f64>>();
    let dim = mat_a.size[0];
    let n = dim as i32;
    let mut w = vec![0.0; dim];
    let mut work = vec![0.0; 4*dim];
    let lwork = 4*n;
    let mut info = 0;


    unsafe {
        dsyev(jobz.clone() as u8, 
              b'L', 
              n, 
              &mut a,
              n, 
              &mut w,
              &mut work,
              lwork, 
              &mut info);
    }

    if info<0 {
        panic!("Error in _dsyev: the {}-th argument had an illegal value", -info);
    } else if info>0 {
        panic!("Error in _dsyev: the algorithm failed to converge; {}-off-diagonal elements of an intermediate tridiagonal form did not converge to zero", info);
    }

    let evec = if jobz == 'V' || jobz == 'v' {
        Some(MatFull::from_vec([dim, dim], a))
    } else {
        None
    };

    (evec, w, n)

}


#[cfg(test)]
mod test {
    use crate::matrix::MatFull;

    use super::{mat_dgemm, mat_dsyev};
    use rand::Rng;

    #[test]
    fn test_dgemm() {
        let mat_a = MatFull::<f64>::from_vec([3,2], vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
        let mat_b = MatFull::<f64>::from_vec([2,3], vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);

        let mat_c = mat_dgemm(&mat_a, &mat_b, 'n','n', 1.0, 1.0);
        assert_eq!(mat_c.data, vec![9.0, 12.0, 15.0, 19.0, 26.0, 33.0, 29.0, 40.0, 51.0]);

        let mat_d = MatFull::<f64>::from_vec([3,3], vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]);
        let mat_e = mat_dgemm(&mat_a, &mat_d, 't','n', 1.0, 1.0);
        assert_eq!(mat_e.data, vec![14.0, 32.0, 32.0, 77.0, 50.0, 122.0]);
    
    }

    #[test]
    fn test_dsyev() {
        let mat_a = MatFull::<f64>::from_vec([3,3], vec![1.0, 2.0, 3.0, 2.0, 4.0, 5.0, 3.0, 5.0, 6.0]);
        let (ev, w, n) = mat_dsyev(&mat_a, 'V');
        ev.unwrap().formated_output_f64('n', 3, "");
        println!("{:?}", w);
        println!("{:?}", n);


    }

    #[test]
    fn dgemm_bench() {
        let mut rng = rand::thread_rng();
        let mat_a = MatFull::<f64>::from_vec([500,500], vec![rng.gen(); 250000]);
        let mat_b = MatFull::<f64>::from_vec([500,500], vec![rng.gen(); 250000]);
        // timer
        let start = std::time::Instant::now();
        let _c = mat_dgemm(&mat_a, &mat_b, 'n','n', 1.0, 1.0);
        println!("Time elapsed in dgemm: {:?}", start.elapsed());
    }


}