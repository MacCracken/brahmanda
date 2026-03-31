use criterion::{Criterion, black_box, criterion_group, criterion_main};

use brahmanda::{halo, morphology, power_spectrum};
use brahmanda::constants::*;

fn bench_nfw_density(c: &mut Criterion) {
    c.bench_function("nfw_density", |b| {
        b.iter(|| halo::nfw_density(black_box(10.0), black_box(1e7), black_box(20.0)).unwrap())
    });
}

fn bench_sersic_profile(c: &mut Criterion) {
    c.bench_function("sersic_profile", |b| {
        b.iter(|| morphology::sersic_profile(black_box(5.0), black_box(10.0), black_box(4.0)).unwrap())
    });
}

fn bench_growth_factor(c: &mut Criterion) {
    c.bench_function("growth_factor", |b| {
        b.iter(|| power_spectrum::growth_factor(black_box(1.0), black_box(OMEGA_M)).unwrap())
    });
}

fn bench_transfer_function(c: &mut Criterion) {
    c.bench_function("transfer_function_eh", |b| {
        b.iter(|| {
            power_spectrum::transfer_function_eh(
                black_box(0.1),
                black_box(OMEGA_M),
                black_box(OMEGA_B),
                black_box(0.674),
            )
            .unwrap()
        })
    });
}

fn bench_halo_properties(c: &mut Criterion) {
    c.bench_function("HaloProperties::from_mass", |b| {
        b.iter(|| halo::HaloProperties::from_mass(black_box(1e12)).unwrap())
    });
}

criterion_group!(
    benches,
    bench_nfw_density,
    bench_sersic_profile,
    bench_growth_factor,
    bench_transfer_function,
    bench_halo_properties,
);
criterion_main!(benches);
