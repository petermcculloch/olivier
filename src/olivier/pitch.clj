(ns olivier.pitch
  (:require [clojure.math.numeric-tower :as math])
  (:use [olivier.sequences :only (rotations positions rev)]
        [clojure.set :only (intersection difference index select project)]
        [markov.core :as markov]))


;==============================================================================
;                            Pitch Functions
;==============================================================================
;
;  A library of functions for operating on pitches via set theory
;


(defn pitchclass
  "Convert to pitchclass"
  ([x] (pitchclass 12 x))
  ([modn x] (mod x modn)))

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
  ([modn pcs amt] (map (partial pc-transpose modn amt) pcs)))

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
  ([coll n] (let [hs  (apply hash-set (mapv pitchclass coll))
                  pcs (range n)] (mapv #(if (contains? hs %) 1 0) pcs))))

(defn bits->pcs [pc-vec] (filter #(>= % 0) (mapv #(- (* % %2) 1) pc-vec (range 1 (count pc-vec)))))


(defn bit-rotate [n bits] (let [len  (count bits)
                                nmod (mod (+ (mod (- n) len) len) len)]
                            (vec (take len (drop nmod (cycle bits))))))
(defn- neg [x] (- x))

(defn flip "x is a bit" [x] (bit-flip x 0))

(defn bit-complement [bits] (mapv (if #(> 0 %) 0 1) bits))

(defn bit-invert [pc-vec] (vec (rseq pc-vec)))

(defn bit-rotations [pitches subdiv]
  (let [pcs      (mapv (partial pitchclass subdiv) pitches)
        inv-pcs  (mapv pc-invert pitches)
        bits     (pitches->bits pcs subdiv)
        inv-bits (pitches->bits inv-pcs subdiv)
        res      (zipmap
                   (concat (mapv #(bit-rotate (neg %) bits) pcs) (mapv #(bit-rotate % inv-bits) inv-pcs))
                   (concat pcs (map neg inv-pcs)))]
    res))


;==============================================================================
;                       Pitchset Analysis and Operations
;==============================================================================


(defn symmetrical-transpositions [modn pcs]
  (let [pcs            (map (partial pitchclass modn) pcs)
        pcs-set        (apply hash-set pcs)
        transpositions (for [x (range 0 modn)]
                         (mapv #(pitchclass modn (+ x %)) pcs))]
    (filterv #(= (apply hash-set %) pcs-set) transpositions)))

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
   (let [pcs        (sort (map pitchclass pcs))
         intervals  (pcs->intervals modn (map #(- % (apply min pcs)) pcs))
         rots       (concat (rotations intervals) (rotations (reverse intervals)))
         forms      (sort-by (comp last vec) (map intervals->pcs rots))
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
   (let [pcs (map (partial pitchclass modn) pcs)]
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
         rots      (vec (concat (rotations intervals) (rotations (reverse intervals))))]
     (mapv (comp vec (partial intervals->pcs false)) (sort left-packer (vec (distinct (mapv vec rots))))))))

(defn all-transpositions
  "Calculate all transpositions for a pitchset"
  ([pcs] (all-transpositions 12 pcs))
  ([modn pcs]
   (map (fn [x] (sort (mapv #(pitchclass modn (+ % x)) pcs))) (range modn))))


(defn in-common
  "Find common pitches across inversions and transpositions"
  ([pcs-a pcs-b] (in-common 12 pcs-a pcs-b))
  ([modn pcs-a pcs-b]
   (let [set-a  (into #{} pcs-a)
         v      (sort-by (comp - count second)
                         (map
                           #(vector % (intersection set-a (into #{} %)))
                           (rotations-and-inversions modn pcs-b)))
         groups (partition-by (comp count last) v)]
     {:best   {:pc     (ffirst v)
               :shared (second (first v))
               :count  (count (second (first v)))
               }
      :groups groups}
     )))


(defn label-pcs
  "Produces a map which converts pitchclasses into set-numbers.
   Handles edge-cases with chromatic/wholetone/diminished/augmented/tritone collections
   by making sure the lowest pitch-class in the set corresponds to the lowest note in the chord
   when all intervals are equal and sum to 12"
  ([pitches] (label-pcs 12 pitches))
  ([modn pitches]
   (let [pc-orig        (distinct (sort (map pitchclass pitches)))
         pc-transp0     (map #(- % (apply min pc-orig)) pc-orig)
         orig-intervals (vec (pcs->intervals pc-orig))
         pc-normal      (normal-form pc-transp0)
         intervals-norm (vec (pcs->intervals pc-normal))
         rot-norm       (rotations intervals-norm)
         rot-inv        (rotations (rseq intervals-norm))
         is-inverted?   (not-any? #(= orig-intervals (vec %)) rot-norm)
         v              (if is-inverted?
                          (nth (partition (count pc-normal) 1 (cycle (reverse pc-normal)))
                               (first (positions #(= (vec %) orig-intervals) rot-inv)))
                          (nth (partition (count pc-normal) 1 (cycle pc-normal))
                               (first (positions #(= (vec %) orig-intervals) rot-norm))))]
     (apply sorted-map (interleave
                         (if (> (count (symmetrical-transpositions modn pc-normal)) 1)
                           (map (partial pitchclass modn) (reductions + (first pitches) intervals-norm))
                           pc-orig)
                         v)))))


(defn- find-closest [x coll]
  (let [distfn    (fn [y] (map #(math/abs (- y %)) coll))
        top       (apply max coll)
        bottom    (apply min coll)
        distances (sort-by second
                           (map vector coll (distfn x) (map (partial min) (distfn top) (distfn bottom))))
        min-dist  (apply min (map second distances))
        distances (filter #(= min-dist (second %)) distances)
        result    (map first (sort-by (juxt second last) distances))]
    (first result)
    ))

;; TODO Close but not quite there.
(defn lengthen-if-needed [modn a b]
  (if (not= (count a) (count b))
    (let [octaves      (sort (distinct (flatten (map (juxt identity dec)
                                                     (map #(math/floor (/ % modn)) (concat a b))))))
          [a b] (if (> (count a) (count b)) [b a] [a b])    ; rebind so a is the smaller set now
          search-space (apply sorted-set (for [x a oct octaves]
                                           (+ x (* oct modn))))
          missing      (difference (apply hash-set (map #(find-closest % search-space) b)) (apply hash-set a))
          results      (sort (concat a (map #(find-closest % search-space) missing)))
          ]
      results) a))


(defn create-rule
  "Where a and b are chords of equal length"
  ([a b] (create-rule 12 a b))
  ([modn a b]
   (let [a    (lengthen-if-needed modn a b)
         b    (lengthen-if-needed modn b a)
         a-pc (map (label-pcs a) (map pitchclass a))
         b-pc (map (label-pcs b) (map pitchclass b))
         dist (map - b a)
         chan (range (count b))]                            ; Uses zero-based counting, rather than Cope's 1-based
     (mapv #(vector [% %2] %3 %4) a-pc b-pc dist chan))))

(defn synthesize-rule
  "Uses voice-leading from rule-a on pitch material from rule-b"
  [rule-a rule-b]
  (let [sa        (sort-by ffirst rule-a)
        sb        (sort-by ffirst rule-b)
        idxs      (map last (sort-by first (partition 2 (map last (interleave sb sa)))))
        edit-item (fn [coll new-idx] (assoc coll 2 new-idx)) ; 2 is index of last value
        new-rule  (sort-by last (mapv #(edit-item % %2) rule-b idxs))]
    new-rule))

(defn apply-chord-transition
  "Given chord-a and a rule for the voice-leading for the pitchsets of a->b,
   produce chord-b"
  [chord-a rule-b]
  (mapv #(+ % (second %2)) chord-a rule-b))

(defn extract-pitchsets-from-rule
  [rule]
  (vector (mapv ffirst rule)
          (mapv (comp second first) rule)))

(defn- make-chord-keys [chord]
  (let [pcs (map pitchclass chord)
        v   (mapv #((label-pcs pcs) %) pcs)
        sv  (vec (distinct (sort v)))]
    {:sorted sv :voiced v :played chord}))


(defn create-rule-map [rules]
  (apply hash-set (for [r rules]
                    (let [[a-pcs b-pcs] (extract-pitchsets-from-rule r)
                          [sa-pcs sb-pcs] (mapv (comp vec sort) [a-pcs b-pcs])]
                      {:pitchset-start (vec (distinct sa-pcs)),
                       :pitchset-end   (vec (distinct sb-pcs)),
                       :voiced-start   (vec a-pcs),
                       :voiced-end     (vec b-pcs),
                       :rule           r}))))


(defn find-voicing-to-any-chords [rule-set voiced-a]
  (select #(= (:voiced-start %) voiced-a) rule-set))


(defn find-pitchset-combinations-with-voicing [rule-map pcs-a pcs-b voicing-a]
  (get-in rule-map [{:pitchset-start pcs-a} {:pitchset-end pcs-b} {:voiced-start voicing-a}]))


(defn find-pitchset-combinations [rule-map pcs-a pcs-b]
  (get-in rule-map [{:pitchset-start pcs-a} {:pitchset-end pcs-b}]))


(defn find-all-voicings-of-chords [rule-map pcs-a pcs-b]
  (let [chord-combinations (find-pitchset-combinations rule-map pcs-a pcs-b)]
    (mapcat vals (merge (apply concat (map #(for [kvs %] (apply hash-map kvs)) (vals chord-combinations)))))))

(defn find-pitchset-combinations-with-voicing [rule-set pcs-a pcs-b voiced-a]
  (select #(and (= (:voiced-start %) voiced-a)
                (= (:pitchset-start pcs-a))
                (= (:pitchset-end pcs-b))) rule-set))


(defn find-pitchset-combinations [rule-set pcs-a pcs-b]
  (select #(and (= (:pitchset-start pcs-a))
                (= (:pitchset-end pcs-b))) rule-set))


(defn find-all-voicings-of-chords [rule-set pcs-a pcs-b]
  (select #(and (= (:pitchset-start pcs-a))
                (= (:pitchset-end pcs-b)))
          rule-set))


(defn get-rule-for-chords [transition-state]
  (let [{:keys [rules start-chord target-pitchset]} transition-state
        rule-set  (create-rule-map rules)
        {sorted-a :sorted voiced-a :voiced} (make-chord-keys start-chord)
        {sorted-b :sorted} (make-chord-keys target-pitchset)
        ; If there is existing voicing, use that chord
        ; Otherwise build a rule with
        ; 1. the exact voicing of the first chord to the next pitch set OR any chord
        ; 2. any voicing of the first chord to the second chord
        voicings1 (select #(= (:voiced-start %) voiced-a) rule-set)
        voicing1  (if (seq voicings1)
                    (if-let [exact (seq (select #(= (:pitchset-end %) sorted-b) voicings1))]
                      (:rule (first exact)) (:rule (rand-nth (vec voicings1)))) nil)
        voicings2 (select #(and (= (:pitchset-start %) sorted-a)
                                (= (:pitchset-end %) sorted-b)) rule-set)
        voicing2  (if (seq voicings2)
                    (if-let [exact (seq (select #(= (:voiced-start %) voiced-a) voicings2))]
                      (:rule (first exact)) (:rule (rand-nth (vec voicings2)))) nil)
        rule1->2  (if (= voicing1 voicing2) voicing1 (synthesize-rule voicing1 voicing2))
        rule2->1  (if (= voicing1 voicing2) voicing1 (synthesize-rule voicing1 voicing2))
        new-rules (if (and (not= rule1->2 rule2->1) (every? not-empty [rule1->2 rule2->1]))
                    [rule1->2 rule2->1] [])
        ]
    (assoc transition-state :start-chord start-chord,
                            :target-pitchset target-pitchset,
                            :rule-start voicing1,
                            :rule-end voicing2,
                            :new-rules new-rules
                            :rule-next rule1->2
                            :rules (if (not-empty new-rules) (apply conj rules new-rules) rules))))





; ======================================================================================================================
;                                               Testing Code
; ======================================================================================================================

(defn build-chord-transition-database [progs]
  (partition 2 1 progs))

(defn build-chord-rules [chord-transitions]
  (apply sorted-set (map (partial apply create-rule) chord-transitions)))

(def chromatic-progressions
  [[48 55 64 72] [48 57 65 72] [50 57 65 72] [50 55 65 71] [52 59 62 71] [52 57 60 69] [50 57 60 65] [52 56 59 64]
   [52 57 60 69] [54 57 64 69] [55 59 62 67] [52 59 68 71] [57 60 69 72] [57 61 67 73] [53 62 65 74] [56 62 65 74]
   [57 64 67 72] [57 63 67 72] [57 62 66 72] [57 62 65 72] [56 62 65 71] [55 62 65 71] [55 60 64 71] [54 60 64 71]
   [54 60 64 69] [54 59 63 69] [52 59 64 68] [50 59 62 65] [50 59 62 67] [48 55 64 67] [52 60 67 72] [53 60 69 72]
   [53 62 69 74] [52 59 68 76] [57 64 67 72] [50 57 69 77] [52 59 68 76] [53 60 69 72] [50 57 72 81] [52 60 72 79]
   [53 62 69 72] [55 59 67 74] [52 60 67 72] [51 60 67 72] [50 60 66 72] [49 60 65 72] [48 55 64 72] [53 60 64 69]
   [52 59 62 67] [50 57 60 65] [48 55 59 64] [47 53 57 62] [52 56 59 64] [48 57 60 64] [49 57 64 67] [50 57 62 65]
   [50 59 62 65] [48 53 60 64] [47 53 57 62] [52 56 59 62] [57 61 61 64] [59 62 62 66] [57 61 64 68] [59 62 66 69]
   [61 64 68 71] [62 66 69 73] [62 65 69 72] [55 62 65 71] [55 60 64 71] [54 60 64 69] [54 59 63 69] [52 59 62 68]
   [52 58 62 67] [52 57 61 67] [50 57 60 65] [50 56 60 65] [50 55 59 65] [48 55 59 64]])

(def played-progressions
  [[48 55 60 64]
   [50 57 60 65] [52 59 62 67] [53 60 64 69] [52 59 62 67] [52 57 60 67] [57 64 67 72] [59 62 67 74] [60 64 67 76]
   [53 60 69 76] [53 62 69 76] [53 62 69 74] [52 60 67 74] [52 60 67 71] [52 60 67 72] [53 60 69 72] [53 62 69 72]
   [55 62 67 71] [59 62 67 71] [53 62 69 77] [50 57 69 77] [52 60 67 76] [57 64 67 76] [53 62 69 74] [55 62 71 74]
   [57 64 67 72] [52 60 64 67] [53 60 65 69] [50 57 65 72] [55 62 67 71] [53 62 67 71] [52 60 67 76] [53 62 69 77]
   [55 64 67 71] [57 64 67 72] [59 65 69 74] [59 64 68 74] [57 64 67 72] [52 60 67 72] [53 60 69 72] [55 64 67 72]
   [56 64 71 74] [57 64 69 72] [56 64 71 74] [55 60 72 76] [54 64 69 72] [55 60 67 72] [53 62 72 74] [52 60 72 76]
   [53 60 69 72] [55 60 67 72] [55 62 67 71] [53 62 67 71] [52 60 67 71] [53 60 64 69] [52 59 62 67] [50 57 60 65]
   [48 55 59 64] [53 60 64 69] [48 55 59 64] [53 60 64 69] [57 65 69 76] [55 64 67 74] [57 62 65 72] [55 60 67 71]
   [57 64 67 72] [57 62 65 72] [55 62 65 71] [48 60 64 72] [48 55 67 72] [48 57 60 64] [48 53 57 60] [48 52 55 60]
   ])



(def counter (atom 0))                                      ; for stepping through chord sequence

(defn initialize-transition
  ([progressions] (initialize-transition progressions (first progressions)))
  ([progressions start-chord] (initialize-transition progressions
                                                     start-chord
                                                     (progressions (inc (.indexOf progressions start-chord)))))
  ([progressions start-chord next-chord]
   (do (reset! counter 0)
       (let [pc-weights (markov/build-from-coll (mapv (fn [xs] (mapv pitchclass xs)) progressions))
             pc-seq     (lazy-seq (markov/generate-walk pc-weights))]
         {
          :pitchset-seq    pc-seq
          :pitchset-gen    (fn [] (do (swap! counter inc) (nth pc-seq (dec @counter))))
          :start-chord     start-chord,
          :target-pitchset (map pitchclass next-chord),
          :rule-next       nil,
          :rule-start      nil,
          :rule-end        nil,
          :rules           (build-chord-rules (build-chord-transition-database played-progressions))}))))


(defn iterate-chords [transition-state]
  (let [t (get-rule-for-chords transition-state)
        {:keys [start-chord rule-next chordlist pitchset-gen]} t
        next-chord (apply-chord-transition start-chord rule-next)]
    (assoc t :start-chord (if (empty next-chord) start-chord next-chord),
             :target-pitchset (pitchset-gen)                ; next-pitchset)
             :chordlist (conj chordlist next-chord))))


(defn demo-chord-output
  ([] (demo-chord-output chromatic-progressions))
  ([progressions] (demo-chord-output progressions 300))
  ([progressions n]
   (let [m      (initialize-transition progressions)
         chords (:chordlist (nth (iterate iterate-chords m) n))]
     (filterv not-empty chords)
     )))






