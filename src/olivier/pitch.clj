(ns olivier.pitch
  (:require [clojure.math.numeric-tower :as math]
            [clojure.math.combinatorics :as combo]
            [clojure.zip :as zip])
  (:use [olivier.sequences :only (rotations positions rev)]
            [clojure.set :only (intersection difference)]))

;==============================================================================
;                            Pitch Functions
;==============================================================================
;
;  A library of functions for operating on pitches via set theory
;


(defn pitchclass
  ([x] (pitchclass 12 x))
  ([div x] (mod (+ (mod x div) div) div)))

(defn wrap-interval
  "Wraps an interval to its smallest form (e.g. 7 -> 5 for mod 12 pitches)"
  ([interval] (wrap-interval interval 12))
  ([interval modn]
   (let [x (mod (math/abs interval) modn)]
     (if (> x (/ modn 2)) (- modn x) x)
     )))

(defn transpose-to-zero
  "Transposes a set of pitches/pitchclasses so that the lowest note is 0"
  [pcs] (mapv #(- % (apply min pcs)) pcs))

(defn pitchset-invert [pcs] (rev (mapv #(- (apply max pcs) %) pcs)))



;==============================================================================
;                           Pitchclass Operations
;==============================================================================

(defn pitches->pcs
  "Returns a sorted vector of the distinct pitchclasses in the list.  "
  [pitches] (vec (distinct (sort (map pitchclass pitches)))))

(defn pcs->intervals
  ([pcs] (pcs->intervals 12 pcs))
  ([modn pcs] (let [pc-transp (sort (transpose-to-zero pcs))]
                (map - (concat (rest pc-transp) [modn]) pc-transp))))

(defn intervals->pcs
  ([intervals] (intervals->pcs false intervals))
  ([uselast? intervals] (reductions + 0 (if uselast? intervals (butlast intervals)))))

(defn pc-transpose
  ([pc amt] (pc-transpose 12 pc amt))
  ([modn pc amt] (pitchclass modn (+ amt pc))))

(defn pcs-transpose
  ([pcs amt] (pcs-transpose 12 pcs amt))
  ([modn pcs amt] (map (partial pc-transpose modn amt) pcs )))

(defn pc-invert
  "Invert pitchclass"
  ([pc] (pc-invert pc 12))
  ([subdiv pc] (pitchclass subdiv (- subdiv pc))))

(defn pc-complement [modn pcs]
  (vec (difference (apply sorted-set (range modn)) (apply sorted-set pcs))))

(defn pcs-invert
  ([pcs] (pcs-invert true pcs))
  ([in-place? pcs]
   (if in-place?
     (vec (sort (map #(+ (- (apply max pcs) %) (apply min pcs)) pcs)))
     (sort (map #(- (apply max pcs) %) pcs)))))


;==============================================================================
;                           Bit Operations
;==============================================================================

(defn pitches->bits
  ([coll] (pitches->bits coll 12))
  ([coll n] (let [hs (apply hash-set (mapv pitchclass coll))
                  pcs (range n)] (mapv #(if (contains? hs %) 1 0) pcs))))

(defn bits->pcs [pc-vec] (filter #(>= % 0) (mapv #(- (* % %2) 1) pc-vec (range 1 (count pc-vec)))))


(defn bit-rotate [n bits] (let [len (count bits)
                                nmod (mod (+ (mod (- n) len) len) len)]
                            (vec (take len (drop nmod (cycle bits))))))
(defn- neg [x] (- x))

(defn flip "x is a bit" [x] (bit-flip x 0))

(defn bit-complement [bits] (mapv (if #(> 0 %) 0 1) bits))

(defn bit-invert [pc-vec] (vec (rseq pc-vec)))

(defn bit-rotations [pitches subdiv]
  (let [pcs (mapv (partial pitchclass subdiv) pitches)
        inv-pcs (mapv pc-invert pitches)
        bits (pitches->bits pcs subdiv)
        inv-bits (pitches->bits inv-pcs subdiv)
        res      (zipmap
                   (concat (mapv #(bit-rotate (neg %) bits) pcs) (mapv #(bit-rotate % inv-bits) inv-pcs))
                   (concat pcs (map neg inv-pcs)))]
    res))


;==============================================================================
;                       Pitchset Analysis and Operations
;==============================================================================

(defn list-of-intervals
  ([pcs] (list-of-intervals pcs 12))
  ([pcs modn]
   "Returns a sorted list of all the wrapped-intervals"
   (loop [pc-coll pcs coll []]
     (if pc-coll
       (recur (next pc-coll) (concat coll (map #(wrap-interval (- (first pc-coll) %) modn) (rest pc-coll))))
       (sort coll)))))

(defn normal-form
  ([pcs] (normal-form 12 pcs))
  ([modn pcs]
   (let [pcs (sort (map pitchclass pcs))
         intervals (pcs->intervals modn (map #(- % (apply min pcs)) pcs))
         rots (concat (rotations intervals) (rotations (reverse intervals)))
         forms (sort-by (comp last vec) (map intervals->pcs rots))
         candidates (sort-by vec (first (partition-by last forms)))] (vec (first candidates)))))


(defn- left-packer [a b]
  (if (= (last a) (last b))
    (compare (vec (butlast a)) (vec (butlast b)))
    (compare (last b) (last a))))

(defn left-pack
  ([pcs] (left-pack 12 pcs))
  ([modn pcs]
   (->> pcs
        (pcs->intervals modn)
        rotations
        (sort left-packer)
        first
        (intervals->pcs false))))

(defn normal-form
  ([pcs] (normal-form 12 pcs))
  ([modn pcs]
   (let [pcs (map pitchclass pcs)]
     (vec (first (sort [(vec (left-pack pcs)) (vec (left-pack (pcs-invert pcs)))]))))))

(defn interval-vector
  ([pcs] (interval-vector pcs 12))
  ([pcs modn]
   (let [intervals (frequencies (list-of-intervals pcs modn))]
     (map #(get intervals % 0) (range 1 (inc (/ modn 2)))))))

(defn rotations-and-inversions
  ([pcs] (rotations-and-inversions 12 pcs))
  ([modn pcs]
  (let [intervals (pcs->intervals modn pcs)
        rots (vec (concat (rotations intervals) (rotations (reverse intervals))))]
    (mapv (comp vec (partial intervals->pcs false)) (sort left-packer (vec (distinct (mapv vec rots))))))))

(defn all-transpositions
  "Calculate all transpositions for a pitchset"
  ([pcs] (all-transpositions 12 pcs))
  ([modn pcs]
   (map (fn [x] (sort (mapv #(pitchclass (+ % x)) pcs))) (range modn))))


(defn in-common
  "Find common pitches across inversions and transpositions"
  ([pcs-a pcs-b] (in-common 12 pcs-a pcs-b))
  ([modn pcs-a pcs-b]
  (let [set-a (into #{} pcs-a)
        v (sort-by (comp - count second)
                   (map
                     #(vector % (intersection set-a (into #{} %)))
                     (rotations-and-inversions modn pcs-b)))
        groups (partition-by (comp count last) v)]
    {:best {:pc (ffirst v)
            :shared (second (first v))
            :count (count (second (first v)))
            }
     :groups groups}
    )))


(defn label-pcs
  "Produces a map which converts pitchclasses into set-numbers.
   Handles edge-cases with chromatic/wholetone/diminished/augmented/tritone collections
   by making sure the lowest pitch-class in the set corresponds to the lowest note in the chord
   when

  "
  [pitches]
  (let [pc-orig            (distinct (sort (map pitchclass pitches)))
        pc-transp0         (map #(- % (apply min pc-orig)) pc-orig)
        orig-intervals     (vec (pcs->intervals pc-orig))
        pc-normal          (normal-form pc-transp0)
        intervals-norm     (vec (pcs->intervals pc-normal))
        rot-norm           (rotations intervals-norm)
        rot-inv            (rotations (rseq intervals-norm))
        is-inverted?       (not-any? #(= orig-intervals (vec %)) rot-norm)
        v                  (if is-inverted?
                             (nth (partition (count pc-normal) 1 (cycle (reverse pc-normal)))
                                  (first (positions #(= (vec %) orig-intervals) rot-inv)))
                             (nth (partition (count pc-normal) 1 (cycle pc-normal))
                                  (first (positions #(= (vec %) orig-intervals) rot-norm))))]
    (apply sorted-map (interleave pc-orig v))))

;; Examples drawn from David Cope's "The Algorithmic Composer", page 90.
;; (label-pcs [36 50 69 77])
;; (label-pcs [45 61 64 78])
;; (label-pcs [50 58 67 76])
;; (label-pcs [62 65 68 69 74])

(pcs-invert true [2 3 6 9])

(defn- find-closest [x coll]
  (let [distfn (fn [y] (map #(math/abs (- y %)) coll))
        top (apply max coll)
        bottom (apply min coll)
        distances (sort-by second
                           (map vector coll (distfn x) (map (partial min) (distfn top) (distfn bottom))))
        min-dist (apply min (map second distances))
        distances (filter #(= min-dist (second %)) distances)
        result (map first (sort-by (juxt second last) distances))]
    (first result)
    ))

;; Close but not quite there.
(defn lengthen-if-needed [modn a b]
  (if (not= (count a) (count b))
    (let [octaves (sort (distinct (flatten (map (juxt identity dec)
                                                (map #(math/floor (/ % modn)) (concat a b))))))
          [a b] (if (> (count a) (count b)) [b a] [a b]) ; rebind so a is the smaller set now
          search-space (apply sorted-set (for [x a oct octaves]
                                           (+ x (* oct modn))))
          missing (difference (apply hash-set (map #(find-closest % search-space) b)) (apply hash-set a))
          results (sort (concat a (map #(find-closest % search-space) missing)))
          ]
      results) a))

;(defn create-rule
;  "Where a and b are chords of equal length"
;  [a b]
;
;  (let [a (lengthen-if-needed 12 a b)
;        b (lengthen-if-needed 12 a b)
;        a-pc (map (label-pcs a) (map pitchclass a))
;        b-pc (map (label-pcs b) (map pitchclass b))
;        dist (map - b a)
;        chan (range (count b))] ; Uses zero-based counting, rather than Cope's 1-based
;    (mapv #(vector [% %2] %3 %4) a-pc b-pc dist chan)))

(defn create-rule
  "Where a and b are chords of equal length"
  ([a b] (create-rule 12 a b))
  ([modn a b]
  (let [a (lengthen-if-needed modn a b)
        b (lengthen-if-needed modn b a)
        a-pc (map (label-pcs a) (map pitchclass a))
        b-pc (map (label-pcs b) (map pitchclass b))
        dist (map - b a)
        chan (range (count b))] ; Uses zero-based counting, rather than Cope's 1-based
      (mapv #(vector [% %2] %3 %4) a-pc b-pc dist chan))))


(defn- synthesize-rule
  "Uses voice-leading from rule-a on pitch material from rule-b"
  [rule-a rule-b]
  (let [sa (sort-by ffirst rule-a)
        sb (sort-by ffirst rule-b)
        sa-idxs (map last sa)
        idxs (map last (sort-by first (partition 2 (map last (interleave sb sa)))))
        edit-item (fn [coll new-idx] (assoc coll 2 new-idx)) ; 2 is index of last value
        new-rule (sort-by last (mapv #(edit-item % %2) rule-b idxs))]
    new-rule))

(defn- apply-chord-transition
  "Given chord-a and a rule for the voice-leading for the pitchsets of a->b,
   produce chord-b"
  [chord-a rule-b]
  (map #(+ % (second %2)) chord-a rule-b)

  )

(defn extract-pitchsets-from-rule
  [rule]
  (vector (mapv ffirst rule)
          (mapv (comp second first) rule)))

(defn- make-chord-keys [chord]
  (let [pcs (map pitchclass chord)
        v (map #((label-pcs pcs) %) pcs)
        sv (sort v)]
    {:sorted sv :voiced v :played chord}))

(comment
(def c1 [60 63 65 68]) ; p 89
(def c2 [62 65 63 67])
(def c3 [36 50 69 77]) ; p 90
(def c4 [38 48 71 76])
(def c5 [45 61 64 78])
(def c6 [50 58 67 76])
(def c7 [36 50 69 77])
(def c8 [39 48 66 82])
;; return mapping of pitch to pc in normal-vec
(def c9 [29 50 74])

(lengthen-if-needed c8 c9)
(lengthen-if-needed c9 c8)


(create-rule c1 c2)
(create-rule c3 c4)
(create-rule c5 c6)
(create-rule c7 c8))

(def ii7 [50 57 65 72])
(def V7b9 [50 56 65 71])
(def Imaj7 [48 55 64 71])
(def I [48 55 64 72])
(def IVmaj7 [48 53 63 69])
(def iii7 [52 59 62 67])
(def vi7 [52 57 60 67])

(def ii7-alt [53 60 69 74])
(def V7b9-alt [53 59 68 74])

(def rule-ii-V (create-rule ii7 V7b9))
(def rule-ii-V-alt (create-rule ii7-alt V7b9))
(apply-chord-transition ii7-alt (synthesize-rule rule-ii-V-alt rule-ii-V))


(def chord-progressions (partition 2 1 [ii7 V7b9 Imaj7 IVmaj7 iii7 vi7 ii7 Imaj7 ii7 iii7 IVmaj7]))
(def chord-progressions (combo/combinations [ii7 V7b9 Imaj7 IVmaj7 iii7 vi7 ii7-alt I Imaj7 ii7 iii7 IVmaj7] 2))



(def chord-rules (apply sorted-set (map (partial apply create-rule) chord-progressions)))



(defn create-rule-map [rules]
  (let [key-rules (map extract-pitchsets-from-rule rules)]
    (loop [rs rules results {}]
      (if (empty? rs)
        results
        (let [r (first rs)
              [a-pcs b-pcs] (map distinct (extract-pitchsets-from-rule r))
              [sa-pcs sb-pcs] (mapv sort [a-pcs b-pcs])]
          (recur (next rs) (assoc-in results [{:pitchset-start sa-pcs} {:pitchset-end sb-pcs} {:voiced-start a-pcs} {:voiced-end b-pcs}] r)))))))

;
(defn- find-nested
  "From: https://stackoverflow.com/questions/28091305/find-value-of-specific-key-in-nested-map"
  [m k]
  (->> (tree-seq map? vals m)
       (filter map?)
       (some k)))

(defn- find-all-nested
  [m k]
  "From: https://stackoverflow.com/questions/28091305/find-value-of-specific-key-in-nested-map"
  (->> (tree-seq map? vals m)
       (filter map?)
       (keep k)))

(defn map-zipper [m]
  (zip/zipper
    (fn [x] (or (map? x) (map? (nth x 1))))
    (fn [x] (seq (if (map? x) x (nth x 1))))
    (fn [x children]
      (if (map? x)
        (into {} children)
        (assoc x 1 (into {} children))))
    m))

(defn find-voicing-to-any-chords [rule-map voiced-a]
  "Returns a map of rules with :voiced-end as key"
  (loop [z (map-zipper rule-map) coll []]
    (if (zip/end? z)
      (if (empty? coll) coll (apply merge (map second coll)))
      (recur (zip/next z)
             (let [res (get (first (zip/node z)) :voiced-start)]
               (if (or (nil? res) (not= res voiced-a))
                 coll
                 (conj coll (zip/node z))))))))


(defn find-pitchset-combinations [rule-map pcs-a pcs-b]
  (get-in rule-map [{:pitchset-start pcs-a} {:pitchset-end pcs-b}]))

(defn find-all-voicings-of-chords [rule-map pcs-a pcs-b]
  (let [chord-combinations (find-pitchset-combinations rule-map pcs-a pcs-b)]
    (mapcat vals (merge (apply concat (map #(for [kvs %] (apply hash-map kvs)) (vals chord-combinations)))))))

(defn get-rule-for-chords [rules current-chord next-pitch-set  prev-chord]
  (let [rule-map  (create-rule-map rules)
        {sorted-a :sorted voiced-a :voiced} (make-chord-keys current-chord)
        {sorted-b :sorted} (make-chord-keys next-pitch-set)

        ; If there is existing voicing, use that chord
        ; Otherwise build a rule with
        ; 1. the exact voicing of the first chord to the next pitch set OR any chord
        ; 2. any voicing of the first chord to the second chord
        voicings1 (find-voicing-to-any-chords rule-map voiced-a)
        voicing1  (if (empty? voicings1) nil (rand-nth (vals voicings1)))     ; or randomly choose
        voicings2 (find-all-voicings-of-chords rule-map sorted-a sorted-b)
        voicing2  (if (empty? voicings2) nil (rand-nth voicings2))
        ]
    {:rule-start voicing1 :rule-end voicing2 :rule (if (and voicing1 voicing2)
                                                     (synthesize-rule voicing1 voicing2)
                                                     nil)
     }))


(defn- print-wide [m]
  (binding [clojure.pprint/*print-right-margin* 200] (clojure.pprint/pprint m)))

(def chord-rule-map (create-rule-map chord-rules))

(map #((label-pcs ii7) (pitchclass %)) ii7)

