{-# LANGUAGE ForeignFunctionInterface #-}

module MathQ where
import Foreign.C
import Foreign.C.String

foreign import ccall "Absolute_Student_t_Distribution" c_Absolute_Student_t_Distribution :: CDouble -> CInt -> CDouble
foreign import ccall "Absolute_Student_t_Distribution_Large_dof" c_Absolute_Student_t_Distribution_Large_dof :: CDouble -> CInt -> CDouble 
foreign import ccall "Auxiliary_Cos_Integral_gi" c_Auxiliary_Cos_Integral_gi :: CDouble -> CDouble 
foreign import ccall "Auxiliary_Sin_Integral_fi" c_Auxiliary_Sin_Integral_fi :: CDouble -> CDouble 
foreign import ccall "Bernoulli_Number" c_Bernoulli_Number :: CInt -> CDouble 
foreign import ccall "Bernoulli_Random_Variate" c_Bernoulli_Random_Variate :: CDouble -> IO CInt
foreign import ccall "Beta_Density" c_Beta_Density :: CDouble -> CDouble -> CDouble -> CDouble
foreign import ccall "Beta_Distribution" c_Beta_Distribution :: CDouble -> CDouble -> CDouble -> CDouble
foreign import ccall "Beta_Function" c_Beta_Function :: CDouble -> CDouble -> CDouble
foreign import ccall "Beta_Random_Variate" c_Beta_Random_Variate :: CDouble -> CDouble -> IO CDouble
foreign import ccall "Binomial_Coefficient" c_Binomial_Coefficient :: CInt -> CInt -> CDouble
foreign import ccall "Binomial_Cumulative_Distribution" c_Binomial_Cumulative_Distribution :: CInt -> CInt -> CDouble -> CDouble
foreign import ccall "Binomial_Point_Distribution" c_Binomial_Point_Distribution :: CInt -> CInt -> CDouble -> CDouble
foreign import ccall "Binomial_Random_Variate" c_Binomial_Random_Variate :: CInt -> CDouble -> IO CInt
foreign import ccall "Catalan_Beta_Function" c_Catalan_Beta_Function :: CDouble -> CDouble
foreign import ccall "Catalan_Beta_Star_Function" c_Catalan_Beta_Star_Function :: CDouble -> CDouble
foreign import ccall "Cauchy_Density" c_Cauchy_Density :: CDouble -> CDouble
foreign import ccall "Cauchy_Distribution" c_Cauchy_Distribution :: CDouble -> CDouble
foreign import ccall "Cauchy_Random_Variate" c_Cauchy_Random_Variate :: IO CDouble
foreign import ccall "Charlier_Cn" c_Charlier_Cn :: CDouble -> CDouble -> CInt -> CDouble
foreign import ccall "Chebyshev_Shifted_Tn" c_Chebyshev_Shifted_Tn :: CDouble -> CInt -> CDouble
foreign import ccall "Chebyshev_Shifted_Un" c_Chebyshev_Shifted_Un :: CDouble -> CInt -> CDouble
foreign import ccall "Chebyshev_Shifted_Vn" c_Chebyshev_Shifted_Vn :: CDouble -> CInt -> CDouble
foreign import ccall "Chebyshev_Shifted_Wn" c_Chebyshev_Shifted_Wn :: CDouble -> CInt -> CDouble
foreign import ccall "Chebyshev_Tn" c_Chebyshev_Tn :: CDouble -> CInt -> CDouble
foreign import ccall "Chebyshev_Un" c_Chebyshev_Un :: CDouble -> CInt -> CDouble
foreign import ccall "Chebyshev_Vn" c_Chebyshev_Vn :: CDouble -> CInt -> CDouble
foreign import ccall "Chebyshev_Wn" c_Chebyshev_Wn :: CDouble -> CInt -> CDouble
foreign import ccall "Chi_Square_Density" c_Chi_Square_Density :: CDouble -> CInt -> CDouble
foreign import ccall "Chi_Square_Distribution" c_Chi_Square_Distribution :: CDouble -> CInt -> CDouble
foreign import ccall "Chi_Square_Distribution_Large_dof" c_Chi_Square_Distribution_Large_dof :: CDouble -> CInt -> CDouble
foreign import ccall "Complete_Elliptic_Integral_First_Kind" c_Complete_Elliptic_Integral_First_Kind :: CChar -> CDouble -> CDouble
foreign import ccall "Complete_Elliptic_Integral_Second_Kind" c_Complete_Elliptic_Integral_Second_Kind :: CChar -> CDouble -> CDouble
foreign import ccall "Cos_Integral_Ci" c_Cos_Integral_Ci :: CDouble -> CDouble
foreign import ccall "Dawsons_Integral" c_Dawsons_Integral :: CDouble -> CDouble
foreign import ccall "DiGamma_Function" c_DiGamma_Function :: CDouble -> CDouble
foreign import ccall "Dirichlet_Eta_Function" c_Dirichlet_Eta_Function :: CDouble -> CDouble
foreign import ccall "Dirichlet_Eta_Star_Function" c_Dirichlet_Eta_Star_Function :: CDouble -> CDouble
foreign import ccall "Dirichlet_Lambda_Function" c_Dirichlet_Lambda_Function :: CDouble -> CDouble
foreign import ccall "Dirichlet_Lambda_Star_Function" c_Dirichlet_Lambda_Star_Function :: CDouble -> CDouble
foreign import ccall "Double_Factorial" c_Double_Factorial :: CInt -> CDouble
foreign import ccall "Entire_Cos_Integral_Cin" c_Entire_Cos_Integral_Cin :: CDouble -> CDouble
foreign import ccall "Entire_Incomplete_Gamma_Function" c_Entire_Incomplete_Gamma_Function :: CDouble -> CDouble -> CDouble
foreign import ccall "Euler_Number" c_Euler_Number :: CInt -> CDouble
foreign import ccall "Exponential_Density" c_Exponential_Density :: CDouble -> CDouble
foreign import ccall "Exponential_Distribution" c_Exponential_Distribution :: CDouble -> CDouble
foreign import ccall "Exponential_Integral_Alpha_n" c_Exponential_Integral_Alpha_n :: CDouble -> CInt -> CDouble
foreign import ccall "Exponential_Integral_Beta_n" c_Exponential_Integral_Beta_n :: CDouble -> CInt -> CDouble
foreign import ccall "Exponential_Integral_Ei" c_Exponential_Integral_Ei :: CDouble -> CDouble
foreign import ccall "Exponential_Integral_Ein" c_Exponential_Integral_Ein :: CDouble -> CDouble
foreign import ccall "Exponential_Integral_En" c_Exponential_Integral_En :: CDouble -> CInt -> CDouble
foreign import ccall "Exponential_Random_Variate" c_Exponential_Random_Variate :: IO CDouble
foreign import ccall "Exponential_Variate_Inversion" c_Exponential_Variate_Inversion :: IO CDouble
foreign import ccall "Exponential_Variate_Ziggurat" c_Exponential_Variate_Ziggurat :: IO CDouble
foreign import ccall "F_Density" c_F_Density :: CDouble -> CInt -> CInt -> CDouble
foreign import ccall "F_Distribution" c_F_Distribution :: CDouble -> CInt -> CInt -> CDouble
foreign import ccall "F_Distribution_Large_Denominator_dof" c_F_Distribution_Large_Denominator_dof :: CDouble -> CInt -> CInt -> CDouble
foreign import ccall "F_Distribution_Large_dofs" c_F_Distribution_Large_dofs :: CDouble -> CInt -> CInt -> CDouble
foreign import ccall "F_Distribution_Large_Numerator_dof" c_F_Distribution_Large_Numerator_dof :: CDouble -> CInt -> CInt -> CDouble
foreign import ccall "Factorial" c_Factorial :: CInt -> CDouble
foreign import ccall "Fresnel_Auxiliary_Cosine_Integral" c_Fresnel_Auxiliary_Cosine_Integral :: CDouble -> CDouble
foreign import ccall "Fresnel_Auxiliary_Sine_Integral" c_Fresnel_Auxiliary_Sine_Integral :: CDouble -> CDouble
foreign import ccall "Fresnel_Cosine_Integral" c_Fresnel_Cosine_Integral :: CDouble -> CDouble
foreign import ccall "Fresnel_Sine_Integral" c_Fresnel_Sine_Integral :: CDouble -> CDouble
foreign import ccall "Gamma_Density" c_Gamma_Density :: CDouble -> CDouble -> CDouble
foreign import ccall "Gamma_Distribution" c_Gamma_Distribution :: CDouble -> CDouble -> CDouble
foreign import ccall "Gamma_Function" c_Gamma_Function :: CDouble -> CDouble
foreign import ccall "Gamma_Function_Max_Arg" c_Gamma_Function_Max_Arg :: IO CDouble
foreign import ccall "Gamma_Random_Variate" c_Gamma_Random_Variate :: CDouble -> IO CDouble
foreign import ccall "Gaussian_Density" c_Gaussian_Density :: CDouble -> CDouble
foreign import ccall "Gaussian_Distribution" c_Gaussian_Distribution :: CDouble -> CDouble
foreign import ccall "Gaussian_Random_Variate" c_Gaussian_Random_Variate :: IO CDouble
foreign import ccall "Gaussian_Variate_Box_Muller" c_Gaussian_Variate_Box_Muller :: IO CDouble
foreign import ccall "Gaussian_Variate_Marsaglias_Ziggurat" c_Gaussian_Variate_Marsaglias_Ziggurat :: IO CDouble
foreign import ccall "Gaussian_Variate_Polar_Marsaglia" c_Gaussian_Variate_Polar_Marsaglia :: IO CDouble
foreign import ccall "Gaussian_Variate_Sum_12_Uniforms" c_Gaussian_Variate_Sum_12_Uniforms :: IO CDouble
foreign import ccall "Gegenbauer_Cn" c_Gegenbauer_Cn :: CDouble -> CDouble -> CInt -> CDouble
foreign import ccall "Geometric_Cumulative_Distribution" c_Geometric_Cumulative_Distribution :: CInt -> CDouble -> CDouble
foreign import ccall "Geometric_Point_Distribution" c_Geometric_Point_Distribution :: CInt -> CDouble -> CDouble
foreign import ccall "Geometric_Random_Variate" c_Geometric_Random_Variate :: CDouble -> IO CInt
foreign import ccall "Gumbels_Maximum_Density" c_Gumbels_Maximum_Density :: CDouble -> CDouble
foreign import ccall "Gumbels_Maximum_Distribution" c_Gumbels_Maximum_Distribution :: CDouble -> CDouble
foreign import ccall "Gumbels_Maximum_Random_Variate" c_Gumbels_Maximum_Random_Variate :: IO CDouble
foreign import ccall "Gumbels_Minimum_Density" c_Gumbels_Minimum_Density :: CDouble -> CDouble
foreign import ccall "Gumbels_Minimum_Distribution" c_Gumbels_Minimum_Distribution :: CDouble -> CDouble
foreign import ccall "Gumbels_Minimum_Random_Variate" c_Gumbels_Minimum_Random_Variate :: IO CDouble
foreign import ccall "Hermite_Hen" c_Hermite_Hen :: CDouble -> CInt -> CDouble
foreign import ccall "Hermite_Hn" c_Hermite_Hn :: CDouble -> CInt -> CDouble
foreign import ccall "Heumans_Lambda_Naught" c_Heumans_Lambda_Naught :: CDouble -> CDouble -> CDouble
foreign import ccall "Hypergeometric_Cumulative_Distribution" c_Hypergeometric_Cumulative_Distribution :: CInt -> CInt -> CInt -> CInt -> CDouble
foreign import ccall "Hypergeometric_Point_Distribution" c_Hypergeometric_Point_Distribution :: CInt -> CInt -> CInt -> CInt -> CDouble
foreign import ccall "Incomplete_Beta_Function" c_Incomplete_Beta_Function :: CDouble -> CDouble -> CDouble -> CDouble
foreign import ccall "Incomplete_Gamma_Function" c_Incomplete_Gamma_Function :: CDouble -> CDouble -> CDouble
foreign import ccall "Inverse_Jacobi_cn" c_Inverse_Jacobi_cn :: CDouble -> CChar -> CDouble -> CDouble
foreign import ccall "Inverse_Jacobi_dn" c_Inverse_Jacobi_dn :: CDouble -> CChar -> CDouble -> CDouble
foreign import ccall "Inverse_Jacobi_sn" c_Inverse_Jacobi_sn :: CDouble -> CChar -> CDouble -> CDouble
foreign import ccall "Jacobi_am" c_Jacobi_am :: CDouble -> CChar -> CDouble -> CDouble
foreign import ccall "Jacobi_cn" c_Jacobi_cn :: CDouble -> CChar -> CDouble -> CDouble
foreign import ccall "Jacobi_dn" c_Jacobi_dn :: CDouble -> CChar -> CDouble -> CDouble
foreign import ccall "Jacobi_Pn" c_Jacobi_Pn :: CDouble -> CDouble -> CDouble -> CInt -> CDouble
foreign import ccall "Jacobi_sn" c_Jacobi_sn :: CDouble -> CChar -> CDouble -> CDouble
foreign import ccall "Jacobi_Zeta_Function" c_Jacobi_Zeta_Function :: CDouble -> CChar -> CDouble -> CDouble
foreign import ccall "Kolmogorov_Asymptotic_Distribution" c_Kolmogorov_Asymptotic_Distribution :: CDouble -> CInt -> CDouble
foreign import ccall "Krawtchouk_Kn" c_Krawtchouk_Kn :: CDouble -> CDouble -> CInt -> CInt -> CDouble
foreign import ccall "Kumaraswamys_Density" c_Kumaraswamys_Density :: CDouble -> CDouble -> CDouble -> CDouble
foreign import ccall "Kumaraswamys_Distribution" c_Kumaraswamys_Distribution :: CDouble -> CDouble -> CDouble -> CDouble
foreign import ccall "Kumaraswamys_Random_Variate" c_Kumaraswamys_Random_Variate :: CDouble -> CDouble -> IO CDouble
foreign import ccall "Laguerre_Ln" c_Laguerre_Ln :: CDouble -> CInt -> CDouble
foreign import ccall "Laguerre_Ln_alpha" c_Laguerre_Ln_alpha :: CDouble -> CDouble -> CInt -> CDouble
foreign import ccall "Laplace_Density" c_Laplace_Density :: CDouble -> CDouble
foreign import ccall "Laplace_Distribution" c_Laplace_Distribution :: CDouble -> CDouble
foreign import ccall "Laplace_Random_Variate" c_Laplace_Random_Variate :: IO CDouble
foreign import ccall "Legendre_Elliptic_Integral_First_Kind" c_Legendre_Elliptic_Integral_First_Kind :: CDouble -> CChar -> CDouble -> CDouble
foreign import ccall "Legendre_Elliptic_Integral_Second_Kind" c_Legendre_Elliptic_Integral_Second_Kind :: CDouble -> CChar -> CDouble -> CDouble
foreign import ccall "Legendre_Pn" c_Legendre_Pn :: CDouble -> CInt -> CDouble
foreign import ccall "Legendre_Shifted_Pn" c_Legendre_Shifted_Pn :: CDouble -> CInt -> CDouble
foreign import ccall "Ln_Beta_Function" c_Ln_Beta_Function :: CDouble -> CDouble -> CDouble
foreign import ccall "Ln_Factorial" c_Ln_Factorial :: CInt -> CDouble
foreign import ccall "Ln_Gamma_Function" c_Ln_Gamma_Function :: CDouble -> CDouble
foreign import ccall "Log_Series_Cumulative_Distribution" c_Log_Series_Cumulative_Distribution :: CInt -> CDouble -> CDouble
foreign import ccall "Log_Series_Point_Distribution" c_Log_Series_Point_Distribution :: CInt -> CDouble -> CDouble
foreign import ccall "Log_Series_Random_Variate" c_Log_Series_Random_Variate :: CDouble -> IO CInt
foreign import ccall "Logistic_Density" c_Logistic_Density :: CDouble -> CDouble
foreign import ccall "Logistic_Distribution" c_Logistic_Distribution :: CDouble -> CDouble
foreign import ccall "Logistic_Random_Variate" c_Logistic_Random_Variate :: IO CDouble
foreign import ccall "Negative_Binomial_Cumulative_Distribution" c_Negative_Binomial_Cumulative_Distribution :: CInt -> CInt -> CDouble -> CDouble
foreign import ccall "Negative_Binomial_Point_Distribution" c_Negative_Binomial_Point_Distribution :: CInt -> CInt -> CDouble -> CDouble
foreign import ccall "Negative_Binomial_Random_Variate" c_Negative_Binomial_Random_Variate :: CInt -> CDouble -> IO CInt
foreign import ccall "Nome" c_Nome :: CDouble -> CDouble
foreign import ccall "Pareto_Density" c_Pareto_Density :: CDouble -> CDouble -> CDouble
foreign import ccall "Pareto_Distribution" c_Pareto_Distribution :: CDouble -> CDouble -> CDouble
foreign import ccall "Pareto_Random_Variate" c_Pareto_Random_Variate :: CDouble -> IO CDouble
foreign import ccall "Poisson_Cumulative_Distribution" c_Poisson_Cumulative_Distribution :: CInt -> CDouble -> CDouble
foreign import ccall "Poisson_Point_Distribution" c_Poisson_Point_Distribution :: CInt -> CDouble -> CDouble
foreign import ccall "Poisson_Random_Variate" c_Poisson_Random_Variate :: CDouble -> IO CInt
foreign import ccall "Quadruple_Factorial" c_Quadruple_Factorial :: CInt -> CDouble
foreign import ccall "Riemann_Zeta_Function" c_Riemann_Zeta_Function :: CDouble -> CDouble
foreign import ccall "Riemann_Zeta_Star_Function" c_Riemann_Zeta_Star_Function :: CDouble -> CDouble
foreign import ccall "Rising_Factorial" c_Rising_Factorial :: CInt -> CInt -> CDouble
foreign import ccall "Sin_Integral_Si" c_Sin_Integral_Si :: CDouble -> CDouble
foreign import ccall "Student_t_Density" c_Student_t_Density :: CDouble -> CInt -> CDouble
foreign import ccall "Student_t_Distribution" c_Student_t_Distribution :: CDouble -> CInt -> CDouble
foreign import ccall "Student_t_Distribution_Large_dof" c_Student_t_Distribution_Large_dof :: CDouble -> CInt -> CDouble
foreign import ccall "t2_Density" c_t2_Density :: CDouble -> CDouble
foreign import ccall "t2_Distribution" c_t2_Distribution :: CDouble -> CDouble
foreign import ccall "t2_Variate_Inversion" c_t2_Variate_Inversion :: IO CDouble
foreign import ccall "Triple_Factorial" c_Triple_Factorial :: CInt -> CDouble
foreign import ccall "Uniform_0_1_Init_Seed" c_Uniform_0_1_Init_Seed :: CInt -> IO ()
foreign import ccall "Uniform_0_1_Init_Time" c_Uniform_0_1_Init_Time :: IO ()
foreign import ccall "Uniform_0_1_Density" c_Uniform_0_1_Density :: CDouble -> CDouble
foreign import ccall "Uniform_0_1_Distribution" c_Uniform_0_1_Distribution :: CDouble -> CDouble
foreign import ccall "Uniform_0_1_Random_Variate" c_Uniform_0_1_Random_Variate :: IO CDouble
foreign import ccall "Weibull_Density" c_Weibull_Density :: CDouble -> CDouble -> CDouble
foreign import ccall "Weibull_Distribution" c_Weibull_Distribution :: CDouble -> CDouble -> CDouble
foreign import ccall "Weibull_Random_Variate" c_Weibull_Random_Variate :: CDouble -> IO CDouble

absolute_student_t_distribution :: Double -> Int -> Double
absolute_student_t_distribution x n = realToFrac $ c_Absolute_Student_t_Distribution (realToFrac x) (fromIntegral n)

absolute_student_t_distribution_large_dof :: Double -> Int -> Double
absolute_student_t_distribution_large_dof x n = realToFrac $ c_Absolute_Student_t_Distribution_Large_dof (realToFrac x) (fromIntegral n)

auxiliary_cos_integral_gi :: Double -> Double 
auxiliary_cos_integral_gi x = realToFrac $ c_Auxiliary_Cos_Integral_gi (realToFrac x)

auxiliary_sin_integral_fi :: Double -> Double 
auxiliary_sin_integral_fi x = realToFrac $ c_Auxiliary_Sin_Integral_fi (realToFrac x)

bernoulli_number :: Int -> Double 
bernoulli_number n = realToFrac $ c_Bernoulli_Number (fromIntegral n)

bernoulli_random_variate :: Double -> IO Int
bernoulli_random_variate p = do return . fromIntegral =<< c_Bernoulli_Random_Variate (realToFrac p)

beta_density :: Double -> Double -> Double -> Double
beta_density x a b = realToFrac $ c_Beta_Density (realToFrac x) (realToFrac a) (realToFrac b)

beta_distribution :: Double -> Double -> Double -> Double
beta_distribution x a b = realToFrac $ c_Beta_Distribution (realToFrac x) (realToFrac a) (realToFrac b)

beta_function :: Double -> Double -> Double
beta_function a b = realToFrac $ c_Beta_Function (realToFrac a) (realToFrac b)

beta_random_variate :: Double -> Double -> IO Double
beta_random_variate a b = do return . realToFrac =<< c_Beta_Random_Variate (realToFrac a) (realToFrac b)

binomial_coefficient :: Int -> Int -> Double
binomial_coefficient n m = realToFrac $ c_Binomial_Coefficient (fromIntegral n) (fromIntegral m)

binomial_cumulative_distribution :: Int -> Int -> Double -> Double
binomial_cumulative_distribution n k p = realToFrac $ c_Binomial_Cumulative_Distribution (fromIntegral n) (fromIntegral k) (realToFrac p)

binomial_point_distribution :: Int -> Int -> Double -> Double
binomial_point_distribution n k p = realToFrac $ c_Binomial_Point_Distribution (fromIntegral n) (fromIntegral k) (realToFrac p)

binomial_random_variate :: Int -> Double -> IO Int
binomial_random_variate n p = do return . fromIntegral =<< c_Binomial_Random_Variate (fromIntegral n) (realToFrac p)

catalan_beta_function :: Double -> Double
catalan_beta_function s = realToFrac $ c_Catalan_Beta_Function (realToFrac s)

catalan_beta_star_function :: Double -> Double
catalan_beta_star_function s = realToFrac $ c_Catalan_Beta_Star_Function (realToFrac s)

cauchy_density :: Double -> Double
cauchy_density x = realToFrac $ c_Cauchy_Density (realToFrac x)

cauchy_distribution :: Double -> Double
cauchy_distribution x = realToFrac $ c_Cauchy_Distribution (realToFrac x)

cauchy_random_variate :: IO Double
cauchy_random_variate = do return . realToFrac =<< c_Cauchy_Random_Variate

charlier_cn :: Double -> Double -> Int -> Double
charlier_cn x a n = realToFrac $ c_Charlier_Cn (realToFrac x) (realToFrac a) (fromIntegral n)

chebyshev_shifted_tn :: Double -> Int -> Double
chebyshev_shifted_tn x n = realToFrac $ c_Chebyshev_Shifted_Tn (realToFrac x) (fromIntegral n)

chebyshev_shifted_un :: Double -> Int -> Double
chebyshev_shifted_un x n = realToFrac $ c_Chebyshev_Shifted_Un (realToFrac x) (fromIntegral n)

chebyshev_shifted_vn :: Double -> Int -> Double
chebyshev_shifted_vn x n = realToFrac $ c_Chebyshev_Shifted_Vn (realToFrac x) (fromIntegral n)

chebyshev_shifted_wn :: Double -> Int -> Double
chebyshev_shifted_wn x n = realToFrac $ c_Chebyshev_Shifted_Wn (realToFrac x) (fromIntegral n)

chebyshev_tn :: Double -> Int -> Double
chebyshev_tn x n = realToFrac $ c_Chebyshev_Tn (realToFrac x) (fromIntegral n)

chebyshev_un :: Double -> Int -> Double
chebyshev_un x n = realToFrac $ c_Chebyshev_Un (realToFrac x) (fromIntegral n)

chebyshev_vn :: Double -> Int -> Double
chebyshev_vn x n = realToFrac $ c_Chebyshev_Vn (realToFrac x) (fromIntegral n)

chebyshev_wn :: Double -> Int -> Double
chebyshev_wn x n = realToFrac $ c_Chebyshev_Wn (realToFrac x) (fromIntegral n)

chi_square_density :: Double -> Int -> Double
chi_square_density x n = realToFrac $ c_Chi_Square_Density (realToFrac x) (fromIntegral n)

chi_square_distribution :: Double -> Int -> Double
chi_square_distribution x n = realToFrac $ c_Chi_Square_Distribution (realToFrac x) (fromIntegral n)

chi_square_distribution_large_dof :: Double -> Int -> Double
chi_square_distribution_large_dof x n = realToFrac $ c_Chi_Square_Distribution_Large_dof (realToFrac x) (fromIntegral n)

complete_elliptic_integral_first_kind :: Char -> Double -> Double
complete_elliptic_integral_first_kind arg x = realToFrac $ c_Complete_Elliptic_Integral_First_Kind (castCharToCChar arg) (realToFrac x)

complete_elliptic_integral_second_kind :: Char -> Double -> Double
complete_elliptic_integral_second_kind arg x = realToFrac $ c_Complete_Elliptic_Integral_Second_Kind (castCharToCChar arg) (realToFrac x)

cos_integral_ci :: Double -> Double
cos_integral_ci x = realToFrac $ c_Cos_Integral_Ci (realToFrac x)

dawsons_integral :: Double -> Double
dawsons_integral x = realToFrac $ c_Dawsons_Integral (realToFrac x)

digamma_function :: Double -> Double
digamma_function x = realToFrac $ c_DiGamma_Function (realToFrac x)

dirichlet_eta_function :: Double -> Double
dirichlet_eta_function s = realToFrac $ c_Dirichlet_Eta_Function (realToFrac s)

dirichlet_eta_star_function :: Double -> Double
dirichlet_eta_star_function s = realToFrac $ c_Dirichlet_Eta_Star_Function (realToFrac s)

dirichlet_lambda_function :: Double -> Double
dirichlet_lambda_function s = realToFrac $ c_Dirichlet_Lambda_Function (realToFrac s)

dirichlet_lambda_star_function :: Double -> Double
dirichlet_lambda_star_function s = realToFrac $ c_Dirichlet_Lambda_Star_Function (realToFrac s)

double_factorial :: Int -> Double
double_factorial n = realToFrac $ c_Double_Factorial (fromIntegral n)

entire_cos_integral_cin :: Double -> Double
entire_cos_integral_cin x = realToFrac $ c_Entire_Cos_Integral_Cin (realToFrac x)

entire_incomplete_gamma_function :: Double -> Double -> Double
entire_incomplete_gamma_function x nu = realToFrac $ c_Entire_Incomplete_Gamma_Function (realToFrac x) (realToFrac nu)

euler_number :: Int -> Double
euler_number n = realToFrac $ c_Euler_Number (fromIntegral n)

exponential_density :: Double -> Double
exponential_density x = realToFrac $ c_Exponential_Density (realToFrac x)

exponential_distribution :: Double -> Double
exponential_distribution x = realToFrac $ c_Exponential_Distribution (realToFrac x)

exponential_integral_alpha_n :: Double -> Int -> Double
exponential_integral_alpha_n x n = realToFrac $ c_Exponential_Integral_Alpha_n (realToFrac x) (fromIntegral n)

exponential_integral_beta_n :: Double -> Int -> Double
exponential_integral_beta_n x n = realToFrac $ c_Exponential_Integral_Beta_n (realToFrac x) (fromIntegral n)

exponential_integral_ei :: Double -> Double
exponential_integral_ei x = realToFrac $ c_Exponential_Integral_Ei (realToFrac x)

exponential_integral_ein :: Double -> Double
exponential_integral_ein x = realToFrac $ c_Exponential_Integral_Ein (realToFrac x)

exponential_integral_en :: Double -> Int -> Double
exponential_integral_en x n = realToFrac $ c_Exponential_Integral_En (realToFrac x) (fromIntegral n)

exponential_random_variate :: IO Double
exponential_random_variate = do return . realToFrac =<< c_Exponential_Random_Variate

exponential_variate_inversion :: IO Double
exponential_variate_inversion = do return . realToFrac =<< c_Exponential_Variate_Inversion

exponential_variate_ziggurat :: IO Double
exponential_variate_ziggurat = do return . realToFrac =<< c_Exponential_Variate_Ziggurat

f_density :: Double -> Int -> Int -> Double
f_density x v1 v2 = realToFrac $ c_F_Density (realToFrac x) (fromIntegral v1) (fromIntegral v2)

f_distribution :: Double -> Int -> Int -> Double
f_distribution f v1 v2 = realToFrac $ c_F_Distribution (realToFrac f) (fromIntegral v1) (fromIntegral v2)

f_distribution_large_denominator_dof :: Double -> Int -> Int -> Double
f_distribution_large_denominator_dof f v1 v2 = realToFrac $ c_F_Distribution_Large_Denominator_dof (realToFrac f) (fromIntegral v1) (fromIntegral v2)

f_distribution_large_dofs :: Double -> Int -> Int -> Double
f_distribution_large_dofs f v1 v2 = realToFrac $ c_F_Distribution_Large_dofs (realToFrac f) (fromIntegral v1) (fromIntegral v2)

f_distribution_large_numerator_dof :: Double -> Int -> Int -> Double
f_distribution_large_numerator_dof f v1 v2 = realToFrac $ c_F_Distribution_Large_Numerator_dof (realToFrac f) (fromIntegral v1) (fromIntegral v2)

factorial :: Int -> Double
factorial n = realToFrac $ c_Factorial (fromIntegral n)

fresnel_auxiliary_cosine_integral :: Double -> Double
fresnel_auxiliary_cosine_integral x = realToFrac $ c_Fresnel_Auxiliary_Cosine_Integral (realToFrac x)

fresnel_auxiliary_sine_integral :: Double -> Double
fresnel_auxiliary_sine_integral x = realToFrac $ c_Fresnel_Auxiliary_Sine_Integral (realToFrac x)

fresnel_cosine_integral :: Double -> Double
fresnel_cosine_integral x = realToFrac $ c_Fresnel_Cosine_Integral (realToFrac x)

fresnel_sine_integral :: Double -> Double
fresnel_sine_integral x = realToFrac $ c_Fresnel_Sine_Integral (realToFrac x)

gamma_density :: Double -> Double -> Double
gamma_density x nu = realToFrac $ c_Gamma_Density (realToFrac x) (realToFrac nu)

gamma_distribution :: Double -> Double -> Double
gamma_distribution x nu = realToFrac $ c_Gamma_Distribution (realToFrac x) (realToFrac nu)

gamma_function :: Double -> Double
gamma_function x = realToFrac $ c_Gamma_Function (realToFrac x)

gamma_function_max_arg :: IO Double
gamma_function_max_arg = do return . realToFrac =<< c_Gamma_Function_Max_Arg

gamma_random_variate :: Double -> IO Double
gamma_random_variate a = do return . realToFrac =<< c_Gamma_Random_Variate (realToFrac a)

gaussian_density :: Double -> Double
gaussian_density x = realToFrac $ c_Gaussian_Density (realToFrac x)

gaussian_distribution :: Double -> Double
gaussian_distribution x = realToFrac $ c_Gaussian_Distribution (realToFrac x)

gaussian_random_variate :: IO Double
gaussian_random_variate = do return . realToFrac =<< c_Gaussian_Random_Variate

gaussian_variate_box_muller :: IO Double
gaussian_variate_box_muller = do return . realToFrac =<< c_Gaussian_Variate_Box_Muller

gaussian_variate_marsaglias_ziggurat :: IO Double
gaussian_variate_marsaglias_ziggurat = do return . realToFrac =<< c_Gaussian_Variate_Marsaglias_Ziggurat

gaussian_variate_polar_marsaglia :: IO Double
gaussian_variate_polar_marsaglia = do return . realToFrac =<< c_Gaussian_Variate_Polar_Marsaglia

gaussian_variate_sum_12_uniforms :: IO Double
gaussian_variate_sum_12_uniforms = do return . realToFrac =<< c_Gaussian_Variate_Sum_12_Uniforms

gegenbauer_cn :: Double -> Double -> Int -> Double
gegenbauer_cn x alpha n = realToFrac $ c_Gegenbauer_Cn (realToFrac x) (realToFrac alpha) (fromIntegral n)

geometric_cumulative_distribution :: Int -> Double -> Double
geometric_cumulative_distribution k p = realToFrac $ c_Geometric_Cumulative_Distribution (fromIntegral k) (realToFrac p)

geometric_point_distribution :: Int -> Double -> Double
geometric_point_distribution k p = realToFrac $ c_Geometric_Point_Distribution (fromIntegral k) (realToFrac p)

geometric_random_variate :: Double -> IO Int
geometric_random_variate p = do return . fromIntegral =<< c_Geometric_Random_Variate (realToFrac p)

gumbels_maximum_density :: Double -> Double
gumbels_maximum_density x = realToFrac $ c_Gumbels_Maximum_Density (realToFrac x)

gumbels_maximum_distribution :: Double -> Double
gumbels_maximum_distribution x = realToFrac $ c_Gumbels_Maximum_Distribution (realToFrac x)

gumbels_maximum_random_variate :: IO Double
gumbels_maximum_random_variate = do return . realToFrac =<< c_Gumbels_Maximum_Random_Variate

gumbels_minimum_density :: Double -> Double
gumbels_minimum_density x = realToFrac $ c_Gumbels_Minimum_Density (realToFrac x)

gumbels_minimum_distribution :: Double -> Double
gumbels_minimum_distribution x = realToFrac $ c_Gumbels_Minimum_Distribution (realToFrac x)

gumbels_minimum_random_variate :: IO Double
gumbels_minimum_random_variate = do return . realToFrac =<< c_Gumbels_Minimum_Random_Variate

hermite_hen :: Double -> Int -> Double
hermite_hen x n = realToFrac $ c_Hermite_Hen (realToFrac x) (fromIntegral n)

hermite_hn :: Double -> Int -> Double
hermite_hn x n = realToFrac $ c_Hermite_Hn (realToFrac x) (fromIntegral n)

heumans_lambda_naught :: Double -> Double -> Double
heumans_lambda_naught amplitude modular_angle = realToFrac $ c_Heumans_Lambda_Naught (realToFrac amplitude) (realToFrac modular_angle)

hypergeometric_cumulative_distribution :: Int -> Int -> Int -> Int -> Double
hypergeometric_cumulative_distribution n1 n2 n k = realToFrac $ c_Hypergeometric_Cumulative_Distribution (fromIntegral n1) (fromIntegral n2) (fromIntegral n) (fromIntegral k)

hypergeometric_point_distribution :: Int -> Int -> Int -> Int -> Double
hypergeometric_point_distribution n1 n2 n k = realToFrac $ c_Hypergeometric_Point_Distribution (fromIntegral n1) (fromIntegral n2) (fromIntegral n) (fromIntegral k)

incomplete_beta_function :: Double -> Double -> Double -> Double
incomplete_beta_function x a b = realToFrac $ c_Incomplete_Beta_Function (realToFrac x) (realToFrac a) (realToFrac b)

incomplete_gamma_function :: Double -> Double -> Double
incomplete_gamma_function x nu = realToFrac $ c_Incomplete_Gamma_Function (realToFrac x) (realToFrac nu)

inverse_jacobi_cn :: Double -> Char -> Double -> Double
inverse_jacobi_cn x arg param = realToFrac $ c_Inverse_Jacobi_cn (realToFrac x) (castCharToCChar arg) (realToFrac param)

inverse_jacobi_dn :: Double -> Char -> Double -> Double
inverse_jacobi_dn x arg param = realToFrac $ c_Inverse_Jacobi_dn (realToFrac x) (castCharToCChar arg) (realToFrac param)

inverse_jacobi_sn :: Double -> Char -> Double -> Double
inverse_jacobi_sn x arg param = realToFrac $ c_Inverse_Jacobi_sn (realToFrac x) (castCharToCChar arg) (realToFrac param)

jacobi_am :: Double -> Char -> Double -> Double
jacobi_am u arg x = realToFrac $ c_Jacobi_am (realToFrac u) (castCharToCChar arg) (realToFrac x)

jacobi_cn :: Double -> Char -> Double -> Double
jacobi_cn u arg x = realToFrac $ c_Jacobi_cn (realToFrac u) (castCharToCChar arg) (realToFrac x)

jacobi_dn :: Double -> Char -> Double -> Double
jacobi_dn u arg x = realToFrac $ c_Jacobi_dn (realToFrac u) (castCharToCChar arg) (realToFrac x)

jacobi_pn :: Double -> Double -> Double -> Int -> Double
jacobi_pn x alpha beta n = realToFrac $ c_Jacobi_Pn (realToFrac x) (realToFrac alpha) (realToFrac beta) (fromIntegral n)

jacobi_sn :: Double -> Char -> Double -> Double
jacobi_sn u arg x = realToFrac $ c_Jacobi_sn (realToFrac u) (castCharToCChar arg) (realToFrac x)

jacobi_zeta_function :: Double -> Char -> Double -> Double
jacobi_zeta_function amplitude arg x = realToFrac $ c_Jacobi_Zeta_Function (realToFrac amplitude) (castCharToCChar arg) (realToFrac x)

kolmogorov_asymptotic_distribution :: Double -> Int -> Double
kolmogorov_asymptotic_distribution dn sample_size = realToFrac $ c_Kolmogorov_Asymptotic_Distribution (realToFrac dn) (fromIntegral sample_size)

krawtchouk_kn :: Double -> Double -> Int -> Int -> Double
krawtchouk_kn x p n k = realToFrac $ c_Krawtchouk_Kn (realToFrac x) (realToFrac p) (fromIntegral n) (fromIntegral k)

kumaraswamys_density :: Double -> Double -> Double -> Double
kumaraswamys_density x a b = realToFrac $ c_Kumaraswamys_Density (realToFrac x) (realToFrac a) (realToFrac b)

kumaraswamys_distribution :: Double -> Double -> Double -> Double
kumaraswamys_distribution x a b = realToFrac $ c_Kumaraswamys_Distribution (realToFrac x) (realToFrac a) (realToFrac b)

kumaraswamys_random_variate :: Double -> Double -> IO Double
kumaraswamys_random_variate a b = do return . realToFrac =<< c_Kumaraswamys_Random_Variate (realToFrac a) (realToFrac b)

laguerre_ln :: Double -> Int -> Double
laguerre_ln x n = realToFrac $ c_Laguerre_Ln (realToFrac x) (fromIntegral n)

laguerre_ln_alpha :: Double -> Double -> Int -> Double
laguerre_ln_alpha x alpha n = realToFrac $ c_Laguerre_Ln_alpha (realToFrac x) (realToFrac alpha) (fromIntegral n)

laplace_density :: Double -> Double
laplace_density x = realToFrac $ c_Laplace_Density (realToFrac x)

laplace_distribution :: Double -> Double
laplace_distribution x = realToFrac $ c_Laplace_Distribution (realToFrac x)

laplace_random_variate :: IO Double
laplace_random_variate = do return . realToFrac =<< c_Laplace_Random_Variate

legendre_elliptic_integral_first_kind :: Double -> Char -> Double -> Double
legendre_elliptic_integral_first_kind amplitude arg x = realToFrac $ c_Legendre_Elliptic_Integral_First_Kind (realToFrac amplitude) (castCharToCChar arg) (realToFrac x)

legendre_elliptic_integral_second_kind :: Double -> Char -> Double -> Double
legendre_elliptic_integral_second_kind amplitude arg x = realToFrac $ c_Legendre_Elliptic_Integral_Second_Kind (realToFrac amplitude) (castCharToCChar arg) (realToFrac x)

legendre_pn :: Double -> Int -> Double
legendre_pn x n = realToFrac $ c_Legendre_Pn (realToFrac x) (fromIntegral n)

legendre_shifted_pn :: Double -> Int -> Double
legendre_shifted_pn x n = realToFrac $ c_Legendre_Shifted_Pn (realToFrac x) (fromIntegral n)

ln_beta_function :: Double -> Double -> Double
ln_beta_function a b = realToFrac $ c_Ln_Beta_Function (realToFrac a) (realToFrac b)

ln_factorial :: Int -> Double
ln_factorial n = realToFrac $ c_Ln_Factorial (fromIntegral n)

ln_gamma_function :: Double -> Double
ln_gamma_function x = realToFrac $ c_Ln_Gamma_Function (realToFrac x)

log_series_cumulative_distribution :: Int -> Double -> Double
log_series_cumulative_distribution k p = realToFrac $ c_Log_Series_Cumulative_Distribution (fromIntegral k) (realToFrac p)

log_series_point_distribution :: Int -> Double -> Double
log_series_point_distribution k p = realToFrac $ c_Log_Series_Point_Distribution (fromIntegral k) (realToFrac p)

log_series_random_variate :: Double -> IO Int
log_series_random_variate p = do return . fromIntegral =<< c_Log_Series_Random_Variate (realToFrac p)

logistic_density :: Double -> Double
logistic_density x = realToFrac $ c_Logistic_Density (realToFrac x)

logistic_distribution :: Double -> Double
logistic_distribution x = realToFrac $ c_Logistic_Distribution (realToFrac x)

logistic_random_variate :: IO Double
logistic_random_variate = do return . realToFrac =<< c_Logistic_Random_Variate

negative_binomial_cumulative_distribution :: Int -> Int -> Double -> Double
negative_binomial_cumulative_distribution n k p = realToFrac $ c_Negative_Binomial_Cumulative_Distribution (fromIntegral n) (fromIntegral k) (realToFrac p)

negative_binomial_point_distribution :: Int -> Int -> Double -> Double
negative_binomial_point_distribution n k p = realToFrac $ c_Negative_Binomial_Point_Distribution (fromIntegral n) (fromIntegral k) (realToFrac p)

negative_binomial_random_variate :: Int -> Double -> IO Int
negative_binomial_random_variate n p = do return . fromIntegral =<< c_Negative_Binomial_Random_Variate (fromIntegral n) (realToFrac p)

nome :: Double -> Double
nome k = realToFrac $ c_Nome (realToFrac k)

pareto_density :: Double -> Double -> Double
pareto_density x a = realToFrac $ c_Pareto_Density (realToFrac x) (realToFrac a)

pareto_distribution :: Double -> Double -> Double
pareto_distribution x a = realToFrac $ c_Pareto_Distribution (realToFrac x) (realToFrac a)

pareto_random_variate :: Double -> IO Double
pareto_random_variate a = do return . realToFrac =<< c_Pareto_Random_Variate (realToFrac a)

poisson_cumulative_distribution :: Int -> Double -> Double
poisson_cumulative_distribution k mu = realToFrac $ c_Poisson_Cumulative_Distribution (fromIntegral k) (realToFrac mu)

poisson_point_distribution :: Int -> Double -> Double
poisson_point_distribution k mu = realToFrac $ c_Poisson_Point_Distribution (fromIntegral k) (realToFrac mu)

poisson_random_variate :: Double -> IO Int
poisson_random_variate mu = do return . fromIntegral =<< c_Poisson_Random_Variate (realToFrac mu)

quadruple_factorial :: Int -> Double
quadruple_factorial n = realToFrac $ c_Quadruple_Factorial (fromIntegral n)

riemann_zeta_function :: Double -> Double
riemann_zeta_function s = realToFrac $ c_Riemann_Zeta_Function (realToFrac s)

riemann_zeta_star_function :: Double -> Double
riemann_zeta_star_function s = realToFrac $ c_Riemann_Zeta_Star_Function (realToFrac s)

rising_factorial :: Int -> Int -> Double
rising_factorial n m = realToFrac $ c_Rising_Factorial (fromIntegral n) (fromIntegral m)

sin_integral_si :: Double -> Double
sin_integral_si x = realToFrac $ c_Sin_Integral_Si (realToFrac x)

student_t_density :: Double -> Int -> Double
student_t_density x n = realToFrac $ c_Student_t_Density (realToFrac x) (fromIntegral n)

student_t_distribution :: Double -> Int -> Double
student_t_distribution x n = realToFrac $ c_Student_t_Distribution (realToFrac x) (fromIntegral n)

student_t_distribution_large_dof :: Double -> Int -> Double
student_t_distribution_large_dof x n = realToFrac $ c_Student_t_Distribution_Large_dof (realToFrac x) (fromIntegral n)

t2_density :: Double -> Double
t2_density x = realToFrac $ c_t2_Density (realToFrac x)

t2_distribution :: Double -> Double
t2_distribution x = realToFrac $ c_t2_Distribution (realToFrac x)

t2_variate_inversion :: IO Double
t2_variate_inversion = do return . realToFrac =<< c_t2_Variate_Inversion

triple_factorial :: Int -> Double
triple_factorial n = realToFrac $ c_Triple_Factorial (fromIntegral n)

uniform_0_1_init_seed :: Int -> IO ()
uniform_0_1_init_seed seed = c_Uniform_0_1_Init_Seed $ fromIntegral seed

uniform_0_1_init_time :: IO ()
uniform_0_1_init_time = c_Uniform_0_1_Init_Time

uniform_0_1_density :: Double -> Double
uniform_0_1_density x = realToFrac $ c_Uniform_0_1_Density (realToFrac x)

uniform_0_1_distribution :: Double -> Double
uniform_0_1_distribution x = realToFrac $ c_Uniform_0_1_Distribution (realToFrac x)

uniform_0_1_random_variate :: IO Double
uniform_0_1_random_variate = do return . realToFrac =<< c_Uniform_0_1_Random_Variate

weibull_density :: Double -> Double -> Double
weibull_density x a = realToFrac $ c_Weibull_Density (realToFrac x) (realToFrac a)

weibull_distribution :: Double -> Double -> Double
weibull_distribution x a = realToFrac $ c_Weibull_Distribution (realToFrac x) (realToFrac a)

weibull_random_variate :: Double -> IO Double
weibull_random_variate a = do return . realToFrac =<< c_Weibull_Random_Variate (realToFrac a)
