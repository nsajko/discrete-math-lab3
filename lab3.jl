# Copyright 2021 Neven Sajko. All rights reserved.
#
# Lab3.prChrNum(f) calculates and prints the chromatic number of
# the simple graph given in the file with the path f.
#
# Ref. https://thorehusfeldt.files.wordpress.com/2010/08/gca.pdf

module Lab3

export prChrNum

# A graph/network represented as an adjacency list.
const Gr = Vector{BitSet}

# sub is interpreted as a bit vector containing the vertices which induce
# a subgraph of g. This checks that the subgraph is
# an independent vertex set
function independentSet(g::Gr, sub::Int)::Bool
	local i::Int = 1
	local s::Int = sub

	while count_ones(s) != 1
		(s&1 != 0) && for adj in g[i]
			((sub >>> (adj - 1))&1 != 0) && (return false)
		end

		i += 1
		s >>>= 1
	end

	true
end

# Calculates the chromatic number using potentially great amounts of time and space.
#
# The algorithm is the dynamic programming algorithm based on iterating through all
# non-empty independent vertex sets of vertex-induced subgraphs.
#
# Currently graphs with more than about 30 vertices are unsupported.
function chromaticNumberSlow(g::Gr)::Int
	local n::Int = size(g, 1)

	(32 < n) && (return 1000) ### XXX TODO: some other algorithm?

	# The table which will contain the chromatic numbers of all
	# vertex-induced subgraphs.
	local cns::Vector{Int8} = zeros(Int8, 1 << n)

	# Mark all one-vertex subgraphs as having chromatic number 1.
	for i in 1 : (1 << (n - 1))
		(count_ones(i) == 1) && (cns[i + 1] = 1)
	end

	for m in 2:n
		for i in 0 : ((1 << n) - 1)
			(count_ones(i) != m) && continue

			local min::Int = 127

			for j in 1:i
				((i|j != i) || !independentSet(g, j)) && continue

				local cn::Int = cns[(i & ~j) + 1]
				(cn < min) && (min = cn)
			end

			cns[i + 1] = 1 + min
		end
	end

	cns[1 << n]
end

function isComplete(g::Gr)::Bool
	local n::Int = size(g, 1) - 1

	for v in g
		(length(v) != n) && (return false)
	end

	true
end

function isCycle(g::Gr)::Bool
	for v in g
		(length(v) != 2) && (return false)
	end

	true
end

function chromaticNumberComponents(s::Vector{Gr})::Int
	# TODO: The facts that each component's chromatic number is independent from
	# the others' and that we only need the maximum could be used to prevent the
	# calculation of some chromatic numbers in some cases. For example, if
	# Brooks' theorem tells us that a component has a chromatic number no
	# greater than 20 and we already know that another component has a chromatic
	# number of 21, then we need not calculate the chromatic number of
	# the former component.

	# Number of components.
	local m::Int = size(s, 1)

	# Maximal chromatic number from among the components.
	local max::Int = 1

	for g in s
		local cn::Int = 0

		# Number of vertices in the graph.
		local n::Int = size(g, 1)

		if isCycle(g)
			# The component is a cycle, thus its chromatic number follows
			# directly from its number of vertices.
			cn = n&1 + 2
		elseif isComplete(g)
			# The graph is complete, thus its chromatic number equals
			# its number of vertices.
			cn = n
		else
			# TODO: check for bipartiteness first (can be done in O(n)).
			cn = chromaticNumberSlow(g)
		end

		max < cn && (max = cn)
	end

	max
end

function dfs(g::Gr, vert::Int, visited::BitSet)::Nothing
	push!(visited, vert)

	for v in setdiff(g[vert], visited)
		dfs(g, v, visited)
	end

	return nothing
end

# Returns the component of the graph g that consists of vertices in verts.
function extractedComponentGraph(g::Gr, verts::BitSet)::Gr
	local n::Int = size(g, 1)
	local m::Int = length(verts)
	local ret::Gr = BitSet[BitSet() for _ in 1:m]
	local oldToNew::Vector{Int16} = zeros(Int16, n)
	local cnt::Int = 0

	for v in verts
		if oldToNew[v] == 0
			cnt = cnt + 1
			oldToNew[v] = cnt
		end

		for adj in g[v]
			if oldToNew[adj] == 0
				cnt = cnt + 1
				oldToNew[adj] = cnt
			end

			push!(ret[oldToNew[v]], oldToNew[adj])
		end
	end

	ret
end

# Takes a graph, returns its nonempty components as separate graphs.
function components(g::Gr)::Vector{Gr}
	local n::Int = size(g, 1)
	local ret::Vector{Gr} = Gr[]
	local visited::BitSet = BitSet()

	for v in 1:n
		if v in visited || isempty(g[v])
			continue
		end
		push!(visited, v)

		local componentVerts::BitSet = BitSet()
		dfs(g, v, componentVerts)
		push!(ret, extractedComponentGraph(g, componentVerts))
	end

	ret
end

# Returns the minimal number of colors necessary to color the graph,
# i.e. the chromatic number.
function chromaticNumber(g::Gr)::Int
	size(g, 1) == 0 && return 0
	return chromaticNumberComponents(components(g))
end

# Warning: no error handling
function readInput(fn::String)::Gr
	local fs::IOStream = open(fn, lock = false)

	local n::Int = parse(Int, readline(fs))
	readline(fs)

	local g::Gr = BitSet[BitSet() for _ in 1:n]
	for i in 1:n
		local line::String = readline(fs)
		for j in 1:n
			(Int(line[2*j - 1]) & 0xf) != 0 && push!(g[i], j)
		end
	end

	g
end

function prChrNum(file::String)::Nothing
	println(chromaticNumber(readInput(file)))

	return nothing
end

end
