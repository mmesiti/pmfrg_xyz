using SpinFRGLattices.SquareLattice

function split_sitesum_v2_noblocking(siteSum, Npairs, Nsum)
    unrolled_kj_m_xk = Vector{@NamedTuple{kj::Int8, m::Int8, xk::Int8}}(undef, sum(Nsum))

    counter = 1
    for Rij = 1:Npairs, ki = 1:Npairs, x in (@view siteSum[:, Rij])
        if x.ki == ki
            unrolled_kj_m_xk[counter] = (kj = x.kj, m = x.m, xk = x.xk)
            counter += 1
        end
    end

    ns = Vector{Int8}(undef, (Npairs * Npairs))

    counter = 1
    for Rij = 1:Npairs, ki = 1:Npairs
        ns[counter] = count(x -> x.ki == ki, (@view siteSum[:, Rij]))
        counter += 1
    end


    return unrolled_kj_m_xk, cumsum(ns)


end



geom = getSquareLattice(16)
unrolled_kj_m_xk, cumsum_ns =
    split_sitesum_v2_noblocking(geom.siteSum, geom.Npairs, geom.Nsum)

start(Rij, ki) =
    let i = ki + (Rij - 1) * geom.Npairs
        if i > 1
            cumsum_ns[i-1] + 1
        else
            1
        end
    end
stop(Rij, ki) = cumsum_ns[ki+(Rij-1)*geom.Npairs]
blockstart(iblock) = round(Int8, blocksize * (iblock - 1)) + 1
blockstop(iblock) = round(Int8, min(blocksize * iblock, geom.Npairs))

nblocks = cld(geom.Npairs, 8) # max_blocksize=8
blocksize = geom.Npairs / nblocks

CacheEntry = @NamedTuple{Type::Int8, Index::Int8}

function search_and_rotate_cache!(cache::Vector{CacheEntry}, what::CacheEntry)
    if what in cache  # cache hit!
        idx = findnext(x -> x == what, cache, 1)
        # rotate top-idx part of cache
        cache[2:idx] = cache[1:idx-1]
        cache[1] = what
        return true
    else # cache miss
        # rotate whole cache, drop last
        cache[2:end] = cache[1:end-1]
        cache[1] = what
        return false
    end
end

"""
Model cache misses of non-blocking approach:
Loop over all values of Rij,
and inner loop over all values of ki.
"""
function model_cache_misses_noblocking(size::Int)
    cache = Vector{CacheEntry}(undef, size)
    n_misses_Rij, n_misses_ki, n_misses_kj = 0, 0, 0
    for Rij = 1:geom.Npairs, ki = 1:geom.Npairs
        for k_spl = start(Rij, ki):stop(Rij, ki)
            (; kj) = unrolled_kj_m_xk[k_spl]
            if !search_and_rotate_cache!(cache, (Type = Int8(1), Index = Int8(Rij)))
                n_misses_Rij += 1
            end
            if !search_and_rotate_cache!(cache, (Type = Int8(2), Index = Int8(ki)))
                n_misses_ki += 1
            end
            if !search_and_rotate_cache!(cache, (Type = Int8(3), Index = Int8(kj)))
                n_misses_kj += 1
            end
        end
    end
    n_misses_Rij, n_misses_ki, n_misses_kj
end

function model_cache_misses(size::Int)
    cache = Vector{CacheEntry}(undef, size)
    n_misses_Rij, n_misses_ki, n_misses_kj = 0, 0, 0
    for Rijblock = 1:nblocks,
        kiblock = 1:nblocks,
        Rij = blockstart(Rijblock):blockstop(Rijblock),
        ki = blockstart(kiblock):blockstop(kiblock)

        for k_spl = start(Rij, ki):stop(Rij, ki)
            (; kj) = unrolled_kj_m_xk[k_spl]
            if !search_and_rotate_cache!(cache, (Type = Int8(1), Index = Int8(Rij)))
                n_misses_Rij += 1
            end
            if !search_and_rotate_cache!(cache, (Type = Int8(2), Index = Int8(ki)))
                n_misses_ki += 1
            end
            if !search_and_rotate_cache!(cache, (Type = Int8(3), Index = Int8(kj)))
                n_misses_kj += 1
            end
        end

    end
    n_misses_Rij, n_misses_ki, n_misses_kj
end

function model_cache_misses_kj_blocking(size::Int)
    cache = Vector{CacheEntry}(undef, size)
    n_misses_Rij, n_misses_ki, n_misses_kj = 0, 0, 0
    for Rijblock = 1:nblocks,
        kiblock = 1:nblocks,
        kjblock = 1:nblocks,
        Rij = blockstart(Rijblock):blockstop(Rijblock),
        ki = blockstart(kiblock):blockstop(kiblock)

        kjlow = blockstart(kjblock)
        kjhigh = blockstop(kjblock)

        for k_spl = start(Rij, ki):stop(Rij, ki)
            (; kj) = unrolled_kj_m_xk[k_spl]
            if kjlow <= kj <= kjhigh
                if !search_and_rotate_cache!(cache, (Type = Int8(1), Index = Int8(Rij)))
                    n_misses_Rij += 1
                end
                if !search_and_rotate_cache!(cache, (Type = Int8(2), Index = Int8(ki)))
                    n_misses_ki += 1
                end
                if !search_and_rotate_cache!(cache, (Type = Int8(3), Index = Int8(kj)))
                    n_misses_kj += 1
                end
            end
        end

    end
    n_misses_Rij, n_misses_ki, n_misses_kj
end

function split_sitesum_v3_allblocked(siteSum, Npairs, Nsum)

    nblocks = cld(Npairs, 8) # max_blocksize=8
    blocksize = Npairs / nblocks
    blockstart(iblock) = round(Int8, blocksize * (iblock - 1)) + 1
    blockstop(iblock) = round(Int8, min(blocksize * iblock, Npairs))

    unrolled_Rij_ki_kj_m_xk =
        Vector{@NamedTuple{Rij::Int8, ki::Int8, kj::Int8, m::Int8, xk::Int8}}(
            undef,
            sum(Nsum),
        )

    counter = 1
    for Rijblock = 1:nblocks,
        kiblock = 1:nblocks,
        kjblock = 1:nblocks,
        Rij = blockstart(Rijblock):blockstop(Rijblock)

        for x in (@view siteSum[:, Rij])
            if x.ki in blockstart(kiblock):blockstop(kiblock) &&
               x.kj in blockstart(kjblock):blockstop(kjblock)

                unrolled_Rij_ki_kj_m_xk[counter] =
                    (Rij = Rij, ki = x.ki, kj = x.kj, m = x.m, xk = x.xk)
                counter += 1
            end
        end
    end

    return unrolled_Rij_ki_kj_m_xk

end

function model_cache_misses_allblocked(size::Int, geom)
    unrolled_Rij_ki_kj_m_xk =
        split_sitesum_v3_allblocked(geom.siteSum, geom.Npairs, geom.Nsum)
    cache = Vector{CacheEntry}(undef, size)
    n_misses_Rij, n_misses_ki, n_misses_kj = 0, 0, 0

    for (; Rij, ki, kj) in unrolled_Rij_ki_kj_m_xk
        if !search_and_rotate_cache!(cache, (Type = Int8(1), Index = Int8(Rij)))
            n_misses_Rij += 1
        end
        if !search_and_rotate_cache!(cache, (Type = Int8(2), Index = Int8(ki)))
            n_misses_ki += 1
        end
        if !search_and_rotate_cache!(cache, (Type = Int8(3), Index = Int8(kj)))
            n_misses_kj += 1
        end
    end
    n_misses_Rij, n_misses_ki, n_misses_kj
end




for size = 18:2:22
    head = "cache size (total \"rows\"): $size"
    println("="^length(head))
    println(head)
    println("no blocking scheme:")
    misses = model_cache_misses_noblocking(size)
    println(misses)
    println("sum:", sum(misses))

    println("blocking scheme:")
    misses = model_cache_misses(size)
    println(misses)
    println("sum:", sum(misses))

    println("KJ blocking scheme:")
    misses = model_cache_misses_kj_blocking(size)
    println(misses)
    println("sum:", sum(misses))

    println("ALL blocking scheme:")
    misses = model_cache_misses_allblocked(size, geom)
    println(misses)
    println("sum:", sum(misses))


    println("="^length(head))

end
