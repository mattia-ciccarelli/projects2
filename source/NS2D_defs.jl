
mutable struct ST6
  mat::Int64
  nodes::Vector{Int64}
end

mutable struct SB3
  nodes::Vector{Int64}
  loadxy::Vector{Float64}
  press::Float64
end

# define a new type
mutable struct Snode
  coor::Vector{Float64}
  udof::Vector{Int64}
  u::Vector{Float64}
  pdof::Int64
  p::Float64
end

mutable struct Sinfo
  mesh_file::String
  NN::Int64
  NE::Int64
  NL::Int64
  uneq::Int64
  pneq::Int64
  neq::Int64
  nK::Int64
  nA::Int64
  #
  nmm::Int64    # number of master modes
  Lmm::Vector{Int64}  # list of mm
  Lmmcj::Vector{Int64}  # list of mm cpnj
  LambdaSym::Vector{Float64}  # list of mm
  nza::Int64    # autonomous
  nzna::Int64    # autonomous
  nrom::Int64    # autonomous
  nMat::Int64    # autonomous
  Ffreq::Int64
  #Fmodes::Vector{Int64}
  #Fmult::Vector{Float64}
  alpha::Float64
  beta::Float64
  neig::Int64   # number of computed modes
  style::Char
  max_order::Int64
  max_orderNA::Int64
  tol::Float64
  Sinfo() = new()
end

