pub mod one2many {
    pub mod prover;
    pub mod verifier;
}
pub mod avss {
    pub mod dealer;
    pub mod party;
}

#[cfg(test)]
mod tests {
    use crate::avss::dealer::Dealer;
    use crate::avss::party::AvssParty;
    use std::mem::size_of;
    use util::algebra::field::Field;
    use util::algebra::polynomial::MultilinearPolynomial;
    use util::{merkle_tree::MERKLE_ROOT_SIZE, random_oracle::RandomOracle};

    use util::algebra::coset::Coset;
    use util::algebra::field::mersenne61_ext::Mersenne61Ext;
    use util::split_n;

    use util::{CODE_RATE, SECURITY_BITS};

    fn output_proof_size(log_n: usize, terminate_round: usize) -> usize {
        let log_t = log_n - 1;
        let oracle = RandomOracle::new(log_t - terminate_round, SECURITY_BITS / CODE_RATE);
        let mut interpolate_cosets = vec![Coset::new(
            1 << (log_t + CODE_RATE),
            Mersenne61Ext::random_element(),
        )];
        for i in 1..log_t {
            interpolate_cosets.push(interpolate_cosets[i - 1].pow(2));
        }
        let polynomial = MultilinearPolynomial::random_polynomial(log_t);

        let x_shift = Mersenne61Ext::random_element();
        let coset_x = Coset::new(1 << log_n, x_shift);
        let mut folding_parameter = vec![];
        let v = split_n((1 << log_t) - 1);
        for i in &v {
            folding_parameter.push(coset_x.pow(*i).all_elements());
        }
        let mut parties = vec![];
        for i in 0..(1 << log_n) {
            let mut open_point = vec![];
            for j in 0..log_t {
                open_point.push(folding_parameter[j][i % folding_parameter[j].len()]);
            }
            parties.push(AvssParty::new(
                log_t - terminate_round,
                &interpolate_cosets,
                open_point,
                &oracle,
            ));
        }
        let mut dealer = Dealer::new(
            log_t - terminate_round,
            &polynomial,
            &interpolate_cosets,
            &oracle,
            &folding_parameter,
        );
        dealer.send_evaluations(&mut parties);
        dealer.commit_functions(&parties);
        dealer.prove();
        dealer.commit_foldings(&parties);
        let (folding, function) = dealer.query();
        let mut folding0 = vec![];
        let mut function0 = vec![];
        for i in 0..(log_t - terminate_round) {
            if i < log_t - terminate_round - 1 {
                folding0.push(folding[i][0].clone());
            }
            function0.push(function[i][0].clone());
        }
        assert!(parties[0].verify(&folding0, &function0));
        folding0.iter().map(|x| x.proof_size()).sum::<usize>()
            + function0.iter().map(|x| x.proof_size()).sum::<usize>()
            + (log_t - terminate_round) * MERKLE_ROOT_SIZE * 2
            + ((1 << terminate_round) + 1) * size_of::<Mersenne61Ext>() * 2
    }

    #[test]
    fn test_proof_size() {
        for i in 10..=21 {
            let proof_size = output_proof_size(i, 1);
            println!(
                "vss proof size of 2^{} parties is {} bytes",
                i,
                proof_size
            );
        }
    }
}
