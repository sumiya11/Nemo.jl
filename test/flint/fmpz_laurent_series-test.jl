@testset "ZZLaurentSeriesRingElem.constructors" begin
  R, x = laurent_series_ring(ZZ, 30, "x")

  @test elem_type(R) == ZZLaurentSeriesRingElem
  @test elem_type(ZZLaurentSeriesRing) == ZZLaurentSeriesRingElem
  @test parent_type(ZZLaurentSeriesRingElem) == ZZLaurentSeriesRing

  @test isa(R, ZZLaurentSeriesRing)

  a1 = x^3 + 2x + 1

  @test isa(a1, ZZLaurentSeriesRingElem)

  b1 = R(a1)

  @test isa(b1, ZZLaurentSeriesRingElem)

  c1 = R(ZZRingElem[1, 3, 5], 3, 5, 0, 1)

  @test isa(c1, ZZLaurentSeriesRingElem)

  g1 = R(1)
  h1 = R(ZZ(2))
  k1 = R()

  @test isa(g1, ZZLaurentSeriesRingElem)
  @test isa(h1, ZZLaurentSeriesRingElem)
  @test isa(k1, ZZLaurentSeriesRingElem)
end

@testset "ZZLaurentSeriesRingElem.printing" begin
  R, x = laurent_series_ring(ZZ, 30, "x")

  @test !occursin(r"{", string(R))

  @test occursin(r"x", string(x^-1 + 1 - x + x^2 + x^5))
end

@testset "ZZLaurentSeriesRingElem.rand" begin
  R, x = laurent_series_ring(ZZ, 10, "x")

  test_rand(R, -12:12, -10:10)
  test_rand(R, -12:12, make(ZZ, -10:10))
end

@testset "ZZLaurentSeriesRingElem.manipulation" begin
  S, x = laurent_series_ring(ZZ, 30, "x")

  @test max_precision(S) == 30

  a = 2x + x^3
  b = O(x^4)

  @test pol_length(a) == 2
  @test pol_length(b) == 0

  @test valuation(a) == 1
  @test valuation(b) == 4

  @test precision(a) == 31
  @test precision(b) == 4

  @test is_gen(gen(S))

  @test iszero(zero(S))

  @test isone(one(S))

  @test is_unit(-1 + x + 2x^2)

  @test isequal(deepcopy(a), a)
  @test isequal(deepcopy(b), b)

  @test coeff(a, 1) == 2
  @test coeff(b, 7) == 0

  @test characteristic(S) == 0
end

@testset "ZZLaurentSeriesRingElem.unary_ops" begin
  R, x = laurent_series_ring(ZZ, 10, "x")
  for iter = 1:300
    f = rand(R, -12:12, -10:10)

    @test isequal(-(-f), f)
    @test iszero(f + (-f))
  end
end

@testset "ZZLaurentSeriesRingElem.binary_ops" begin
  R, x = laurent_series_ring(ZZ, 10, "x")
  for iter = 1:100
    f = rand(R, -12:12, -10:10)
    g = rand(R, -12:12, -10:10)
    h = rand(R, -12:12, -10:10)
    @test isequal(f + g, g + f)
    @test isequal(f + (g + h), (f + g) + h)
    @test isequal(f*g, g*f)
    @test isequal(f*(g*h), (f*g)*h)
    @test isequal(f - g, -(g - f))
    @test (f - h) + h == f
    @test f*(g + h) == f*g + f*h
    @test f*(g - h) == f*g - f*h
  end
end

@testset "ZZLaurentSeriesRingElem.adhoc_binary_ops" begin
  R, x = laurent_series_ring(ZZ, 10, "x")
  for iter = 1:500
    f = rand(R, -12:12, -10:10)
    c1 = rand(ZZ, -10:10)
    c2 = rand(ZZ, -10:10)
    d1 = rand(zz, -10:10)
    d2 = rand(zz, -10:10)

    @test isequal(c1*f - c2*f, (c1 - c2)*f)
    @test isequal(c1*f + c2*f, (c1 + c2)*f)
    @test isequal(d1*f - d2*f, (d1 - d2)*f)
    @test isequal(d1*f + d2*f, (d1 + d2)*f)

    @test isequal(f*c1 - f*c2, f*(c1 - c2))
    @test isequal(f*c1 + f*c2, f*(c1 + c2))
    @test isequal(f*d1 - f*d2, f*(d1 - d2))
    @test isequal(f*d1 + f*d2, f*(d1 + d2))
  end
end

@testset "ZZLaurentSeriesRingElem.comparison" begin
  R, x = laurent_series_ring(ZZ, 10, "x")
  for iter = 1:500
    f = rand(R, -12:12, -10:10)
    g = deepcopy(f)
    h = R()
    while iszero(h)
      h = rand(R, -12:12, -10:10)
    end

    @test f == g
    @test isequal(f, g)
    @test (precision(h) > min(precision(f), precision(g)) || f != g + h)
    @test (precision(h) > min(precision(f), precision(g)) || !isequal(f, g + h))
  end
end

@testset "ZZLaurentSeriesRingElem.adhoc_comparison" begin
  R, x = laurent_series_ring(ZZ, 10, "x")
  for iter = 1:500
    f = R()
    while f == 0
      f = rand(R, 0:0, -10:10)
    end
    f += rand(R, 1:12, -10:10)
    c1 = rand(ZZ, -10:10)
    d1 = rand(zz, -10:10)

    @test R(c1) == c1
    @test c1 == R(c1)
    @test R(d1) == d1
    @test d1 == R(d1)

    @test R(c1) != c1 + f
    @test c1 != R(c1) + f
    @test R(d1) != d1 + f
    @test d1 != R(d1) + f
  end
end

@testset "ZZLaurentSeriesRingElem.powering" begin
  R, x = laurent_series_ring(ZZ, 10, "x")

  for iter = 1:100
    f = rand(R, -12:12, -10:10)
    r2 = R(1)

    for expn = 0:10
      r1 = f^expn

      @test (f == 0 && expn == 0 && r1 == 0) || isequal(r1, r2)

      r2 *= f
    end
  end
end

@testset "ZZLaurentSeriesRingElem.shift" begin
  R, x = laurent_series_ring(ZZ, 10, "x")
  for iter = 1:300
    f = rand(R, -12:12, -10:10)
    s = rand(0:12)

    @test isequal(shift_right(shift_left(f, s), s), f)
    @test isequal(shift_left(f, s), x^s*f)
    @test precision(shift_right(f, s)) == precision(f) - s
  end
end

@testset "ZZLaurentSeriesRingElem.truncation" begin
  R, x = laurent_series_ring(ZZ, 10, "x")
  for iter = 1:300
    f = rand(R, -12:12, -10:10)
    s = rand(-12:12)

    @test truncate(f, s) == f
    @test isequal(truncate(f, s), f + O(x^s))
    @test precision(truncate(f, s)) == min(precision(f), s)
  end
end

@testset "ZZLaurentSeriesRingElem.inversion" begin
  R, x = laurent_series_ring(ZZ, 10, "x")
  for iter = 1:300
    f = R()
    while iszero(f) || !is_unit(polcoeff(f, 0))
      f = rand(R, -12:12, -10:10)
    end

    @test f*inv(f) == 1
  end
end

@testset "ZZLaurentSeriesRingElem.square_root" begin
  R, x = laurent_series_ring(ZZ, 10, "x")
  for iter = 1:300
    f = rand(R, -12:12, -10:10)
    g = f^2

    @test isequal(sqrt(g)^2, g)
  end
end

@testset "ZZLaurentSeriesRingElem.exact_division" begin
  R, x = laurent_series_ring(ZZ, 10, "x")
  for iter = 1:300
    f = rand(R, -12:12, -10:10)
    g = rand(R, -12:12, -10:10)
    while iszero(g) || !is_unit(polcoeff(g, 0))
      g = rand(R, -12:12, -10:10)
    end

    @test divexact(f, g)*g == f
  end
end

@testset "ZZLaurentSeriesRingElem.adhoc_exact_division" begin
  R, x = laurent_series_ring(ZZ, 10, "x")
  for iter = 1:300
    f = rand(R, -12:12, -10:10)
    c = ZZ()
    while c == 0
      c = rand(ZZ, -10:10)
    end

    @test isequal(divexact(f*c, c), f)
  end
end

@testset "ZZLaurentSeriesRingElem.special_functions" begin
  R, x = laurent_series_ring(ZZ, 10, "x")

  @test isequal(exp(2x + x^2 + O(x^3)), 1 + 2*x + 3*x^2 + O(x^3))
end
