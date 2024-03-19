extern crate criterion;

use criterion::*;

use util::algebra::field::Field;
use util::algebra::polynomial::MultilinearPolynomial;
use util::random_oracle::RandomOracle;
use vss::avss::dealer::Dealer;
use vss::avss::party::AvssParty;

use util::algebra::coset::Coset;
use util::algebra::field::mersenne61_ext::Mersenne61Ext;
use util::split_n;

use util::{CODE_RATE, SECURITY_BITS};
fn vss_deal(log_n: usize, terminate_round: usize) {
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
    dealer.query();
}

fn vss_verify(c: &mut Criterion, log_n: usize, terminate_round: usize) {
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
    let mut group = c.benchmark_group("verify proof");
    group.sample_size(10);
    group.bench_function(format!("vss verify 2^{} parties", log_n), move |b| {
        b.iter(|| {
            parties[0].verify(&folding0, &function0);
        })
    });
}

fn bench_vss_deal(c: &mut Criterion) {
    for i in 10..21 {
        let terminate_round = 1;
        c.bench_function(&format!("vss prove 2^{} parties", i), move |b| {
            b.iter(|| {
                vss_deal(i, terminate_round);
            })
        });
    }
}

fn bench_vss_verify(c: &mut Criterion) {
    for i in 10..21 {
        let terminate_round = 1;
        vss_verify(c, i, terminate_round);
    }
}

criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10);
    targets = bench_vss_deal, bench_vss_verify
);
criterion_main!(benches);
