#include "mlcp.hpp"
/**
 * @file mlcp.cpp
 */
namespace Makena {

using namespace std;


void MLCP::solveSLECholeskyInto_z_1(
    SymMat&                  M,
    const VarVec&            q,
    const long               num,
    vector<enum columnType>& types,
    VarVec&                  z_1
) {
    decompose_cholesky_submatrix(num, types, M);

    solve_lower_diagonal_submatrix(q, num, types, M, z_1);

    solve_upper_diagonal_submatrix(num, types, M, z_1);
}


void MLCP::decompose_cholesky_submatrix(
    const long               num,
    vector<enum columnType>& types,
    SymMat&                  M
) {

    for (long k = 1 ; k <= num; k++){
        auto& tk = types[k-1];
        if ( tk == BILATERAL_FREE              ||
             tk == BILATERAL_BOXED_NOT_CLAMPED ||
             tk == UNILATERAL_ACTIVE              ) {

            M.val(k,k) = sqrt(M.const_val(k,k));

            for(long m = k + 1 ; m <= num ; m++){
                auto& tm = types[m-1];
                if ( tm == BILATERAL_FREE              ||
                     tm == BILATERAL_BOXED_NOT_CLAMPED ||
                     tm == UNILATERAL_ACTIVE              ) {

                    if (fabs(M.const_val(k,k))<EPSILON_LCP_ZERO) {
                        cerr << "M.const_val(" << k << "," << k <<
                                 ")<EPSILON_LCP_ZERO CKP1\n";
                    }
                    M.val(m,k) = M.const_val(m,k) / M.const_val(k,k); 
                }
            }

            for(long j = k + 1; j <= num; j++){
                auto& tj = types[j-1];
                if ( tj == BILATERAL_FREE              ||
                     tj == BILATERAL_BOXED_NOT_CLAMPED ||
                     tj == UNILATERAL_ACTIVE              ) {
                    for(long m = j; m <= num; m++){
                        auto& tm = types[m-1];
                        if ( tm == BILATERAL_FREE              ||
                             tm == BILATERAL_BOXED_NOT_CLAMPED ||
                             tm == UNILATERAL_ACTIVE              ) {
                            M.val(m,j)= M.const_val(m,j) - 
                                        M.const_val(m,k) * M.const_val(j,k) ;
                        }
                    }
                }
            }
        }
    }
}


void MLCP::solve_lower_diagonal_submatrix(
    const VarVec&            q,
    const long               num,
    vector<enum columnType>& types,
    const SymMat&            M,
    VarVec&                  z_1
) {
    // Forward substitution L·y + b = 0 using z as y.
    for(long i = 1; i <= num; i++){
        auto& ti = types[i-1];
        if ( ti == BILATERAL_FREE              ||
             ti == BILATERAL_BOXED_NOT_CLAMPED ||
             ti == UNILATERAL_ACTIVE              ) {
            double sum = 0.0;
            for(long j = 1; j <= i-1; j++){
                auto& tj = types[j-1];
                if ( tj == BILATERAL_FREE              ||
                     tj == BILATERAL_BOXED_NOT_CLAMPED ||
                     tj == UNILATERAL_ACTIVE              ) {
                    sum += (M.const_val(i,j) * z_1.const_val(j));
                }
            }
            if (fabs(M.const_val(i,i))<EPSILON_LCP_ZERO) {
                cerr << "M.const_val(" << i << "," << i << 
                        ")<EPSILON_LCP_ZERO CKP2\n";
            }
            z_1.val(i) = (0.0 - q.const_val(i) - sum)  / M.const_val(i,i);
        }
    }
}


void MLCP::solve_upper_diagonal_submatrix(
    const long               num,
    vector<enum columnType>& types,
    const SymMat&            M,
    VarVec&                  z_1
) {
    // Backward substitution L^t·z - y = 0 using z as y.
    for(long i = num; i >= 1; i--){
        auto& ti = types[i-1];
        if ( ti == BILATERAL_FREE              ||
             ti == BILATERAL_BOXED_NOT_CLAMPED ||
             ti == UNILATERAL_ACTIVE              ) {

            double sum = 0.0;
            for(long j = num; j >= i+1; j--){
                auto& tj = types[j-1];
                if ( tj == BILATERAL_FREE              ||
                     tj == BILATERAL_BOXED_NOT_CLAMPED ||
                     tj == UNILATERAL_ACTIVE              ) {
                    sum += (M.const_val(j,i) * z_1.const_val(j));
                }
            }

            if (fabs(M.const_val(i,i))<EPSILON_LCP_ZERO) {
                cerr << "M.const_val(" << i << "," << i << 
                        ")<EPSILON_LCP_ZERO CKP3\n";
            }
            z_1.val(i) = ( z_1.const_val(i) - sum ) / M.const_val(i,i);
        }
        else {
            z_1.val(i) = 0.0;
        }
    }
}


void MLCP::pgs(
    const SymMat&            M,
    const VarVec&            q,
    VarVec&                  z,
    const VarVec&            z_low,
    const VarVec&            z_high,
    const long               num,
    vector<enum columnType>& types,
    const long               kgs
) {
    for(long r = 0; r < kgs ; r++){

        for(long i = 1; i<= num; i++){
            auto& ti = types[i-1];
            if (ti==BILATERAL_FREE) {
                double delta = q.const_val(i);
                for(long j = 1; j <= num; j++){
                    auto& tj = types[i-1];
                    if (tj != NOT_USED                    && 
                        ti !=BILATERAL_BOXED_CLAMPED_LOW  &&
                        ti !=BILATERAL_BOXED_CLAMPED_HIGH   ) {

                        delta += M.const_val(i,j) * z.const_val(j);
                    }
                }

                if (fabs(M.const_val(i,i))<EPSILON_LCP_ZERO) {
                    cerr << "M.const_val(" << i << "," << i << 
                            ")<EPSILON_LCP_ZERO CKP4\n";
                }
                z.val(i) = z.const_val(i) - delta / M.const_val(i,i);
            }
            else if (ti==BILATERAL_BOXED_NOT_CLAMPED ||
                     ti==BILATERAL_BOXED_CLAMPED_LOW ||
                     ti==BILATERAL_BOXED_CLAMPED_HIGH       ) {

                double delta = q.const_val(i);
                for(long j = 1; j <= num; j++){
                    auto& tj = types[i-1];
                    if (tj != NOT_USED ) {

                        delta += M.const_val(i,j) * z.const_val(j);
                    }
                }

                if (fabs(M.const_val(i,i))<EPSILON_LCP_ZERO) {
                    cerr << "M.const_val(" << i << "," << i << 
                            ")<EPSILON_LCP_ZERO CKP5\n";
                }
                z.val(i) = z.const_val(i) - delta / M.const_val(i,i);
                z.val(i) = std::max(z.const_val(i), z_low.const_val(i));
                z.val(i) = std::min(z.const_val(i), z_high.const_val(i));
            }
            else if (ti == UNILATERAL_NON_ACTIVE||ti == UNILATERAL_ACTIVE) {
                double delta = q.const_val(i);
                for(long j = 1; j <= num; j++){
                    auto& tj = types[i-1];
                    if (tj != NOT_USED                    && 
                        ti !=BILATERAL_BOXED_CLAMPED_LOW  &&
                        ti !=BILATERAL_BOXED_CLAMPED_HIGH   ) {

                        delta += M.const_val(i,j) * z.const_val(j);
                    }
                }

                if (fabs(M.const_val(i,i))<EPSILON_LCP_ZERO) {
                    cerr << "M.const_val(" << i << "," << i << 
                            ")<EPSILON_LCP_ZERO CKP6\n";
                }
                z.val(i) = 

                       std::max(0.0, z.const_val(i) - delta/M.const_val(i,i)); 
            }
        }
    }

    // State update.
    for (long i = 1; i <= num; i++) {
        auto& ti = types[i-1];        
        switch (ti) {

          case BILATERAL_BOXED_NOT_CLAMPED:
          case BILATERAL_BOXED_CLAMPED_LOW:
          case BILATERAL_BOXED_CLAMPED_HIGH:

            if (z.const_val(i) <= z_low.const_val(i)+EPSILON_LCP_ZERO) {
                ti = BILATERAL_BOXED_CLAMPED_LOW;
                z.val(i) = z_low.const_val(i);
            }
            else if (z.const_val(i) >= z_high.const_val(i)-EPSILON_LCP_ZERO) {
                ti = BILATERAL_BOXED_CLAMPED_HIGH;
                z.val(i) = z_high.const_val(i);
            }
            else {
                ti = BILATERAL_BOXED_NOT_CLAMPED;
            }
            break;

          case UNILATERAL_NON_ACTIVE:
          case UNILATERAL_ACTIVE:

            if (z.const_val(i) <= EPSILON_LCP_ZERO) {
                ti = UNILATERAL_NON_ACTIVE;
                z.val(i) = 0.0;
            }
            else {
                ti = UNILATERAL_ACTIVE;
            }
            break;

          default:
            ;
        }
    }
}


long MLCP::solve(
    const SymMat&            M,
    const VarVec&            q,
    VarVec&                  z,
    const VarVec&            z_low,
    const VarVec&            z_high,
    const long               num,
    vector<enum columnType>& types,
    const long               maxIter,
    const long               kgs,
    const long               ksm,
    const double             tolerance,
    WorkMemory&              wm
) {

    double na = 0.0;
    double nb = 0.0;
    long bif=0, bib=0, uni=0;

    for (int i = 1; i <= num; i++) {

        auto& t = types[i-1];
         if ( t==BILATERAL_FREE ) {
            bif++;
            na = std::max(na, fabs(q.const_val(i)));

        }
        else if ( t == BILATERAL_BOXED_NOT_CLAMPED ||
                  t == BILATERAL_BOXED_CLAMPED_LOW ||
                  t == BILATERAL_BOXED_CLAMPED_HIGH  ) {
            bib++;
            if (z.const_val(i) >= z_high.const_val(i)) {
                t = BILATERAL_BOXED_CLAMPED_HIGH;
                z.val(i) = z_high.const_val(i);
            }
            else if (z.const_val(i) <= z_low.const_val(i)) {
                t = BILATERAL_BOXED_CLAMPED_LOW;
                z.val(i) = z_low.const_val(i);
            }
            else {
                t = BILATERAL_BOXED_NOT_CLAMPED;
            }

        }
        else if (t==UNILATERAL_NON_ACTIVE||t==UNILATERAL_ACTIVE) {
            uni++;
            z.val(i) = 0.0;// Reset to zero to avoid pathetic (huge) solution.
            nb = std::max(nb, fabs(q.const_val(i)));

        }
    }

    double coeff_a = 1.0 / (1.0 + na);
    double coeff_b = 1.0 / (1.0 + nb);
    double coeff_c = 1.0 / (1.0 + nb * nb);

    double residue = 0.0;
    for (long i = 0; i < maxIter; i++) {

        pgs(M, q, z, z_low, z_high, num, types, kgs);

        residue = residual(M, q, z, num, types, coeff_a, coeff_b, coeff_c);

        if (residue <= tolerance) {
            return 2*(i+1)-1;
        }

        if(isEligibleForSubspaceMinimization(types, num)) {
            bool feasible = 
                 minimizeSubspace(M, q, z, z_low, z_high, num, types, ksm, wm);
            if (feasible) {
                if (residual(M, q, z, num, types, coeff_a, coeff_b, coeff_c) < 
                                                                   tolerance) {
                    return 2*(i+1);
                }
            }
            else if (i == maxIter-1) {
                // Run pgs one more time before giving up.
                pgs(M, q, z, z_low, z_high, num, types, kgs);
            }
        }
    }
    return maxIter;
}


bool MLCP::isEligibleForSubspaceMinimization(
    vector<enum columnType>& types,
    const long               num
) {
    for (long i = 1; i <= num; i++) {
        auto& ti = types[i-1];    
        if (ti == BILATERAL_BOXED_NOT_CLAMPED || ti == UNILATERAL_ACTIVE) {
            return true;            
        }
    }
    return false;
}


bool MLCP::minimizeSubspace(
    const SymMat&            M,
    const VarVec&            q,
    VarVec&                  z,
    const VarVec&            z_low,
    const VarVec&            z_high,
    const long               num,
    vector<enum columnType>& types,
    const long               kgs,
    WorkMemory&              wm
) {
    long i;
    for (i = 0; i < kgs ; i++) {
        if (i==0) {
            wm.mz_0 = z;
        }
        else if (i==1) {
            wm.mz_0 = wm.mz_b;
        }
        else {
            wm.mz_0 = wm.mz_s;
        }

        wm.mMcopy = M;
        solveSLECholeskyInto_z_1(wm.mMcopy, q, num, types, wm.mz_1);
                               
        bool feasible;

        double alpha = findAlpha_z_0_z_1(
                                       num, types, feasible, wm.mz_0, wm.mz_1);
        if (i==0) {
            if (feasible) {
                z = wm.mz_1;
                // Clamp the values for bilateral boxed constraints.
                for (long i = 1; i <= num; i++) {
                    auto& ti = types[i-1];
                    if (ti==BILATERAL_BOXED_NOT_CLAMPED) {
                        if (z.const_val(i) > z_high.const_val(i)) {
                            feasible = false;
                            ti = BILATERAL_BOXED_CLAMPED_HIGH;
                            z.val(i) = z_high.const_val(i);
                        }
                        else if (z.const_val(i) < z_low.const_val(i)) {
                            feasible = false;
                            ti = BILATERAL_BOXED_CLAMPED_LOW;
                            z.val(i) = z_low.const_val(i);
                        }
                    }
                }
                return feasible;
            }
            interpolate_z_0_z_1(
                  alpha, wm.mz_b, num, wm.mz_0, wm.mz_1, types, z_low, z_high);
        }
        else {
            interpolate_z_0_z_1(
                  alpha, wm.mz_s, num, wm.mz_0, wm.mz_1, types, z_low, z_high);
            if (feasible) {
                break;
            }
        }
    }

    double phi_b = phi(M, q, num, types, wm.mz_b);
    double phi_s = phi(M, q, num, types, wm.mz_s);

    if (phi_b < phi_s) {
        z = wm.mz_b;
    }
    else {
        z = wm.mz_s;
    }

    // Clamp the values for bilateral boxed constraints.
    bool feasible = true;
    for (long i = 1; i <= num; i++) {
        auto& ti = types[i-1];
        if (ti==BILATERAL_BOXED_NOT_CLAMPED) {
            if (z.const_val(i) > z_high.const_val(i)) {
                feasible = false;
                ti = BILATERAL_BOXED_CLAMPED_HIGH;
                z.val(i) = z_high.const_val(i);
            }
            else if (z.const_val(i) < z_low.const_val(i)) {
                feasible = false;
                ti = BILATERAL_BOXED_CLAMPED_LOW;
                z.val(i) = z_low.const_val(i);
            }
        }
        else if (ti==UNILATERAL_NON_ACTIVE) {
            z.val(i) = 0.0;
        }
    }

    return feasible;
}


double MLCP::findAlpha_z_0_z_1(
    const long               num,
    vector<enum columnType>& types,
    bool&                    feasible,
    VarVec&                  z_0,
    VarVec&                  z_1
) {
    feasible = true;
    double alpha  = 1.0;
    for (size_t i = 1; i <= num; i++) {
        auto& ti = types[i-1];
        if (ti == UNILATERAL_ACTIVE) {
            double v1 = z_0.const_val(i);
            double v2 = z_1.const_val(i);
            if (v2 < 0.0 && fabs(v1-v2)>=EPSILON_LCP_ZERO) {
                feasible = false;
                alpha = std::min(alpha, v1/(v1 - v2));
            }
        }
    }
    return alpha;
}


void MLCP::interpolate_z_0_z_1(
    const double             alpha, 
    VarVec&                  zOut,
    const long               num,
    VarVec&                  z_0,
    VarVec&                  z_1,
    vector<enum columnType>& types,
    const VarVec&            z_low,
    const VarVec&            z_high
) {
    double one_minus_alpha = 1.0 - alpha;

    for (size_t i = 1; i <= num; i++) {
        auto& ti = types[i-1];
        switch (ti) {

          case BILATERAL_FREE:
          case BILATERAL_BOXED_NOT_CLAMPED:
            {
                double u1 = z_0.const_val(i);
                double u2 = z_1.const_val(i);
                zOut.val(i) = one_minus_alpha * u1 + alpha * u2;
            }
            break;

          case BILATERAL_BOXED_CLAMPED_LOW:
            zOut.val(i) = z_low.const_val(i);
            break;

          case BILATERAL_BOXED_CLAMPED_HIGH:
            zOut.val(i) = z_high.const_val(i);
            break;

          case UNILATERAL_ACTIVE:
            {
                double u1 = z_0.const_val(i);
                double u2 = z_1.const_val(i);
                // Numerical imprecision may make the value negative.
                double candidate = one_minus_alpha * u1 + alpha * u2;
                if (candidate > EPSILON_LCP_ZERO) {
                    zOut.val(i) = candidate;
                }
                else {
                    // SM has produced negative result. Clamping to zero and
                    // make it NON_ACTIVE. This will be excluded from SM next
                    // time onwards.
                    zOut.val(i) = 0.0;
                    ti =  UNILATERAL_NON_ACTIVE;
                }
            }
            break;

          case UNILATERAL_NON_ACTIVE:
          default:
            zOut.val(i) = 0.0;

        }
    }
}


double MLCP::phi(
    const SymMat&            M,
    const VarVec&            q,
    const long               num,
    vector<enum columnType>& types,
    const VarVec&            z
) {
    double sum_outer = 0.0;

    for(long i = 1; i <= num; i++){
        auto& ti = types[i-1];
        if ( ti == BILATERAL_FREE              ||
             ti == BILATERAL_BOXED_NOT_CLAMPED ||
             ti == UNILATERAL_ACTIVE              ) {

            double sum_inner = 0.0;
            for(long j = 1; j <= num; j++){
                auto& tj = types[j-1];
                if ( tj == BILATERAL_FREE              ||
                     tj == BILATERAL_BOXED_NOT_CLAMPED ||
                     tj == UNILATERAL_ACTIVE              ) {

                    sum_inner += (M.const_val(i,j) * z.const_val(j));

                }
            }

            sum_inner = sum_inner * 0.5;
            sum_inner += q.const_val(i);
            sum_outer += (sum_inner * z.const_val(i));
        }
    }

    return sum_outer;
}


double MLCP::residual(
    const SymMat&            M,
    const VarVec&            q,
    VarVec&                  z,
    const long               num,
    vector<enum columnType>& types,
    const double             coeff_a,
    const double             coeff_b,
    const double             coeff_c
) {
    double max_a = 0.0;
    double max_b = 0.0;
    double max_c = 0.0;

    for(int i = 1; i <= num; i++){
        auto& ti = types[i-1];

        if ( ti == BILATERAL_FREE ) {
            double sum = 0.0;        

            for(int j = 1; j <= num; j++){
                auto& tj = types[j-1];
                if (tj != NOT_USED) {
                    sum += ( M.const_val(i,j) * z.const_val(j) );
                }
            }
            sum += q.const_val(i);
            max_a = std::max(max_a, fabs(sum));
        }
        else if ( ti==UNILATERAL_NON_ACTIVE || ti==UNILATERAL_ACTIVE ) {
            double sum = 0.0;

            for(int j = 1; j <= num; j++){
                auto& tj = types[j-1];
                if (tj != NOT_USED) {
                    sum += ( M.const_val(j,i) * z.const_val(j) ) ;
                }
            }
            sum += q.const_val(i);

            max_b = std::max(max_b, std::min(z.const_val(i), sum));

            max_c = std::max(max_c, std::max(0.0, -1.0*sum));

        }
    }

    return std::max(std::max(fabs(max_a)*coeff_a, 
                             fabs(max_b)*coeff_b),
                             fabs(max_c)*coeff_c);
}


}// namespace Makena
