extern crate criterion;

use criterion::*;

use avss::avss::dealer::Dealer;
use avss::avss::party::AvssParty;
use util::algebra::field::Field;
use util::algebra::polynomial::MultilinearPolynomial;
use util::random_oracle::RandomOracle;

use util::algebra::coset::Coset;
use util::algebra::field::mersenne61_ext::Mersenne61Ext;
use util::split_n;
use util::{CODE_RATE, SECURITY_BITS};

fn avss_deal(log_n: usize, terminate_round: usize) {
    let log_t = log_n - 2;
    let log_d = log_t * 2;
    let oracle = RandomOracle::new(log_d - terminate_round, SECURITY_BITS / CODE_RATE);
    let mut interpolate_cosets = vec![Coset::new(
        1 << (log_t * 2 + CODE_RATE),
        Mersenne61Ext::random_element(),
    )];
    for i in 1..log_d {
        interpolate_cosets.push(interpolate_cosets[i - 1].pow(2));
    }
    let polynomial = MultilinearPolynomial::random_polynomial(log_d);

    let x_shift = Mersenne61Ext::random_element();
    let coset_x = Coset::new(1 << log_n, x_shift);
    let mut folding_parameter = vec![];
    let v = split_n((1 << log_t) - 1);
    for i in &v {
        folding_parameter.push(coset_x.pow(*i).all_elements());
    }
    let y_shift = Mersenne61Ext::random_element();
    let coset_y = Coset::new(1 << log_n, y_shift);
    let last_len = folding_parameter.last().unwrap().len();
    for i in &v {
        folding_parameter.push(
            coset_y
                .pow(*i)
                .all_elements()
                .iter()
                .map(|x| (0..last_len).into_iter().map(|_| *x).collect::<Vec<_>>())
                .flatten()
                .collect(),
        );
    }
    let mut parties = vec![];
    for i in 0..(1 << (log_n * 2)) {
        let mut open_point = vec![];
        for j in 0..log_d {
            open_point.push(folding_parameter[j][i % folding_parameter[j].len()]);
        }
        parties.push(AvssParty::new(
            log_d - terminate_round,
            &interpolate_cosets,
            open_point,
            &oracle,
        ));
    }
    let mut dealer = Dealer::new(
        log_d - terminate_round,
        &polynomial,
        &interpolate_cosets,
        &oracle,
        &folding_parameter,
    );
    dealer.send_evaluations(&mut parties);
    dealer.commit_functions(&parties);
    dealer.prove();
    dealer.commit_foldings(&parties);
    dealer.query();
}

fn bench_avss_deal(c: &mut Criterion) {
    for i in 7..12 {
        let terminate_round = 1;
        c.bench_function(
            &format!(
                "avss prove {} degree and {} parties",
                1 << (i - 2),
                (1 << (i - 2)) * 3 + 1
            ),
            move |b| {
                b.iter(|| {
                    avss_deal(i, terminate_round);
                })
            },
        );
    }
}

fn avss_verify(criterion: &mut Criterion, log_n: usize, terminate_round: usize) {
    let log_t = log_n - 2;
    let log_d = log_t * 2;
    let oracle = RandomOracle::new(log_d - terminate_round, SECURITY_BITS / CODE_RATE);
    let mut interpolate_cosets = vec![Coset::new(
        1 << (log_t * 2 + CODE_RATE),
        Mersenne61Ext::random_element(),
    )];
    for i in 1..log_d {
        interpolate_cosets.push(interpolate_cosets[i - 1].pow(2));
    }
    let polynomial = MultilinearPolynomial::random_polynomial(log_d);

    let x_shift = Mersenne61Ext::random_element();
    let coset_x = Coset::new(1 << log_n, x_shift);
    let mut folding_parameter = vec![];
    let v = split_n((1 << log_t) - 1);
    for i in &v {
        folding_parameter.push(coset_x.pow(*i).all_elements());
    }
    let y_shift = Mersenne61Ext::random_element();
    let coset_y = Coset::new(1 << log_n, y_shift);
    let last_len = folding_parameter.last().unwrap().len();
    for i in &v {
        folding_parameter.push(
            coset_y
                .pow(*i)
                .all_elements()
                .iter()
                .map(|x| (0..last_len).into_iter().map(|_| *x).collect::<Vec<_>>())
                .flatten()
                .collect(),
        );
    }
    let mut parties = vec![];
    for i in 0..(1 << (log_n * 2)) {
        let mut open_point = vec![];
        for j in 0..log_d {
            open_point.push(folding_parameter[j][i % folding_parameter[j].len()]);
        }
        parties.push(AvssParty::new(
            log_d - terminate_round,
            &interpolate_cosets,
            open_point,
            &oracle,
        ));
    }
    let mut dealer = Dealer::new(
        log_d - terminate_round,
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
    for i in 0..(log_d - terminate_round) {
        if i < log_d - terminate_round - 1 {
            folding0.push(folding[i][0].clone());
        }
        function0.push(function[i][0].clone());
    }
    criterion.bench_function(
        &format!(
            "avss verify {} degree and {} parties",
            1 << log_t,
            (1 << log_t) * 3 + 1
        ),
        move |b| {
            b.iter(|| {
                assert!(parties[0].verify(&folding0, &function0));
            })
        },
    );
}

fn bench_avss_verify(c: &mut Criterion) {
    for i in 7..12 {
        let terminate_round = 1;
        avss_verify(c, i, terminate_round);
    }
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = bench_avss_deal, bench_avss_verify
);
criterion_main!(benches);
