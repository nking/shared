package algorithms.sat;

import algorithms.graphs.StronglyConnectedComponents3;

import java.util.*;

public class Boolean2SAT {

    /**
     * given an array of rows of CNF clauses of 2 variables, return whether
     * a solution if the problem is satisfiable, else return null.
     * Each clause is given as the disjunction ("or") of a pair of literals,
     * where a literal can be positive or negative.
     * e.g. a CNF clause {-1, 3} is {@literal (-x1 \/ x3)}.
     The r.t.c. is O(m + n) where m is the number of variables,
     and n is the number of clauses.
     <pre>
     following description in Competetive Programming Handbook by Antti Laaksonen
     Chap 17.2, but adding DPLL pre-processing to solve pure literals and unit variables first.
     </pre>
     * @param cnfClauses an array of rows of CNF clauses of 2 variables.
     * Each clause is given as the disjunction (that is, "or")
     * of a pair of literals,
     * where a literal can be positive or negative.
     * e.g. a CNF clause {-1, 3} is {@literal (-x1 \/ x3)}.
     *                   allowed variable number range is [1, m], inclusive.
     * @param the number of variables.  e.g. for a problem having
     *            {-1,3}, {1,2} there are 3 variables: 1, 2, and 3.
     * @return a map of key = positive literal variable number, value = true or false.
     * note that null is returned for unsatisfiable.
     */
    public Map<Integer, Boolean> solve(int[][] cnfClauses, int m) {

        /*
        (0) Preprocess:
            (A) preprocess for pure literals:
              a pure literal, x_i, appears in clauses as only positive
              or only negative.
              (I) if x_i is a pure literal positive appearing, then
                    it can immed be set to "true",
                  else if x_i is a pure literal negative,
                    it can be immed set to "false".
              (II) simplify the clauses containing the pure literals
                   (have removed this to keep the code simpler to read)
           (B) find clauses containing just 1 variable
               the result of the clause must be "true" so some variables
               can also be set from that.

        (1) create graph:
        for each clause {a,b}, create 2 edges: (-a, b) and (-b, a).

        (2) Find the strongly connected components.

        (3) The problem is solvable if no nodes xi and ¬xi belong to the same
            strongly connected component.

            If a solution exists, the values for the variables can be found by going through the nodes of
            the component graph in a reverse topological sort order.

        (4) create the TS and solve in reverse TS order

        will create a data structure to hold clauses as key=variable in positive form, value=set of variables
        its paired with in positive or negative form as needed.

        x1, -x3
        clause map={(1:-3), (3:-1)}
           as we determine x1=true, we use the map to find x3 = false, and remove both entries after solved.

         */

        // ===== pre-processing ==========
        Map<Integer, Boolean> pureLiterals = findPureLiterals(cnfClauses, m);

        Map<Integer, Boolean> soln = new HashMap<>();
        for (Map.Entry<Integer, Boolean> entry : pureLiterals.entrySet()) {
            soln.put(entry.getKey(), entry.getValue());
        }

        /*
        clauseMap key = literal in positive or negative form and value = paired literal in positive or negative form
         */
        Map<Integer, Set<Integer>> clauseMap = buildClauseMap(soln, cnfClauses, m);
        if (clauseMap == null) {
            return null;
        }

        // ======== create graph from clause Map ====
        //for each clause {a,b}, create 2 edges: (-a, b) and (-b, a).
        Map<Integer, Collection<Integer>> graph = new HashMap<>();
        for (Map.Entry<Integer, Set<Integer>> entry : clauseMap.entrySet()) {
            int a = entry.getKey();
            for (int b : entry.getValue()) {
                graph.putIfAbsent(-a, new HashSet<>());
                graph.get(-a).add(b);

                graph.putIfAbsent(-b, new HashSet<>());
                graph.get(-b).add(a);
            }
        }

        // find the strongly connected components of graph
        StronglyConnectedComponents3 scc = new StronglyConnectedComponents3();
        Map<Integer, Set<Integer>> outputComponents = new HashMap<>();
        Map<Integer, Set<Integer>> sccCondensedAdjMap = new HashMap<>();
        scc.find(graph, outputComponents, sccCondensedAdjMap);

        for (Map.Entry<Integer, Set<Integer>> entry : outputComponents.entrySet()) {
            // if comp contains positive and negative of any literal, it is not solvable
            Set<Integer> comp = entry.getValue();
            for (int c : comp) {
                if (comp.contains(-c)) {
                    return null;
                }
            }
        }

        // solvable, need TS of graph
        List<Integer> ts = topSort(sccCondensedAdjMap);

        if (ts == null) {
            throw new IllegalStateException("Error in algorithm.  clause graph strongly connected components" +
                    " show that problem is solvable, but the topological sort did not include all nodes.");
        }

        for (int i = ts.size() - 1; i >= 0; --i) {
            int root = ts.get(i);
            Set<Integer> comp = outputComponents.get(root);

            for (int u : comp) {

                if (soln.containsKey(Math.abs(u))) continue;

                // solve from existing clauses
                Set<Integer> vs = clauseMap.get(u);

                assert (vs != null);

                for (int v : vs) {
                    if (soln.containsKey(Math.abs(v))) {
                        boolean vV = soln.get(Math.abs(v));
                        // -u /\ -v, -u /\ v, u /\ v, u /\ -v
                        //  only if right side result is F, we have a definite uV to make clause true
                        if ((v < 0 && vV) || (v>0 && !vV)) {
                           if (u < 0) {//(-u, -T) so u is F
                               soln.put(-u, false);
                           } else {
                               soln.put(u, true);
                           }
                           break;
                        }
                    }
                }

                if (soln.containsKey(Math.abs(u))) continue;

                // else set signu to signTrue
                if (u < 0) {
                    soln.put(-u, false);
                } else {
                    soln.put(u, true);
                }
            }// end loop over 1 strongly connected component
        }// end loop over all strongly connected components
        assert(soln.size() == m);
        return soln;
    }

    private List<Integer> topSort(Map<Integer, Set<Integer>> graph) {
        Map<Integer, Integer> incoming = new HashMap<>();
        Set<Integer> vertices = new HashSet<>();
        for (Map.Entry<Integer, Set<Integer>> entry : graph.entrySet()) {
            int u = entry.getKey();
            for (int v : entry.getValue()) {
                incoming.put(v, incoming.getOrDefault(u, 0) + 1);
                vertices.add(v);
            }
            vertices.add(u);
        }
        // vertices with no incoming go in queue
        Queue<Integer> q = new ArrayDeque();
        for (int u : vertices) {
            if (!incoming.containsKey(u)) {
                q.offer(u);
            }
        }

        List<Integer> out = new ArrayList<>();
        while (!q.isEmpty()) {
            int u = q.poll();
            out.add(u);
            if (!graph.containsKey(u)) continue;
            for (int v : graph.get(u)) {
                if (incoming.get(v) == 1) {
                    q.offer(v);
                    incoming.remove(v);
                } else {
                    incoming.put(v, incoming.get(v) - 1);
                }
            }
        }
        if (out.size() < vertices.size()) {
            return null;
        }
        return out;
    }

    // build a clause map
    private Map<Integer, Set<Integer>> buildClauseMap(Map<Integer, Boolean> soln, int[][] cnfClauses, int m) {

        Map<Integer, Set<Integer>> clauseMap = new HashMap<>();

        int a, b;
        for (int[] clause : cnfClauses) {
            a = clause[0];
            b = clause[1];

            /*
            Boolean solvedClause = processUnit(a, b, soln);
            if (solvedClause == null) {
                return null;
            }
            if (solvedClause) continue;
            solvedClause = processUnit(b, a, soln);
            if (solvedClause == null) {
                return null;
            }
            if (solvedClause) continue;
            */
            //add to clause map
            clauseMap.putIfAbsent(a, new HashSet<>());
            clauseMap.get(a).add(b);
            clauseMap.putIfAbsent(b, new HashSet<>());
            clauseMap.get(b).add(a);
        }
        return clauseMap;
    }

    /**
     *
     * @param a
     * @param b
     * @param soln
     * @return True if processed and solved both variables (b added to soln),
     * else if False did not find  a in soln so no effect, else NULL if found a and b solved
     * but they conflict (both are false).
     */
    private Boolean processUnit(int a, int b, Map<Integer, Boolean> soln) {
        if (a < 0) {
            a *= -1;
            b *= -1;
        }
        if (soln.containsKey(a)) {
            boolean aV = soln.get(a);
            if (soln.containsKey(Math.abs(b))) {
                boolean bV = soln.get(Math.abs(b));
                if (b < 0) {
                    bV = !bV;
                }
                if (!aV && !bV) {
                    // not solvable.  conflict of pure literals or unit vars
                    return null;
                }
            } else {
                //(a || -b)  (a || b)  (-a || -b) (-a || b)
                // T   T,F       T,F    F    *F    *F     T
                if (!aV) {// can solve b
                    if (b < 0) {
                        // b is false
                        soln.put(Math.abs(b), false);
                    } else {
                        // b is true
                        soln.put(Math.abs(b), false);
                    }
                    return Boolean.TRUE;
                }
            }
        }
        return Boolean.FALSE;
    }

    private Map<Integer, Boolean> findPureLiterals(int[][] cnfClauses, int m) {
        Set<Integer> posLiterals = new HashSet<>();
        Set<Integer> negLiterals = new HashSet<>();
        for (int[] clause : cnfClauses) {
            if (clause[0] < 0) {
                negLiterals.add(-clause[0]);
            } else {
                posLiterals.add(clause[0]);
            }
            if (clause[1] < 0) {
                negLiterals.add(-clause[1]);
            } else {
                posLiterals.add(clause[1]);
            }
        }
        Map<Integer, Boolean> pureLiterals = new HashMap<>();
        for (int i = 1; i <= m; ++i) {
            if (posLiterals.contains(i) && !negLiterals.contains(i)) {
                pureLiterals.put(i, true);
            } else if (!posLiterals.contains(i) && negLiterals.contains(i)) {
                pureLiterals.put(i, false);
            }
        }
        return pureLiterals;
    }

    /*
    for more than 3 variables in a clause, there are several algorithms including DPLL,
    but none are efficient - the problem is NP-Hard
     */
}
