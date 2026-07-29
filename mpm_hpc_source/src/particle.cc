#include "particle.h"
#include "factory.h"
#include "particle_base.h"
#include "twophase_particle.h"
#include "threephase_particle.h"
#include "particle_threephase_new.h"
#include "particle_threephase_lag.h"
#include "Frozen_particle.h"
#include "Frozen_particle_unsaturated.h"
#include "MHBS_particle.h"

// Single phase particle2D (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::Particle<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2d("P2D");

// Single phase particle3D (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::Particle<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3d("P3D");

// Two phase particle2D (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::TwoPhaseParticle<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2d2phase("P2D2PHASE");

// Two phase particle3D (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::TwoPhaseParticle<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3d2phase("P3D2PHASE");

// Three phase particle2D (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::ThreePhaseParticle<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2d3phase("P2D3PHASE");

// Three phase particle3D (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::ThreePhaseParticle<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3d3phase("P3D3PHASE");

// Three phase particle2D (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::ThreePhaseParticleNew<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2d3phasenew("P2D3PHASE&NEW");

// Three phase particle3D (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::ThreePhaseParticleNew<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3d3phasenew("P3D3PHASE&NEW");

// Three phase particle2D (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::ThreePhaseParticleLag<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2d3phaselag("P2D3PHASE&LAG");

// Three phase particle3D (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::ThreePhaseParticleLag<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3d3phaselag("P3D3PHASE&LAG");

// Two phase particle2D with phase change (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::FrozenParticle<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2d2phasePC("P2D2PHASE&PC");

// Two phase particle3D with phase change (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::FrozenParticle<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3d2phasePC("P3D2PHASE&PC");

// Three phase particle2D with phase change (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::UnsatFrozenParticle<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2d3phasePC("P2D3PHASE&PC");

// Three phase particle3D with phase change (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::UnsatFrozenParticle<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3d3phasePC("P3D3PHASE&PC");

// Two phase particle2D with phase change (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::MHBSParticle<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    MHBSparticle2d("P2DMHBS");

// Two phase particle3D with phase change (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::MHBSParticle<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    MHBSparticle3d("P3DMHBS");