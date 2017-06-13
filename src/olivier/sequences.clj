(ns olivier.sequences
  (:require [clojure.math.combinatorics :as combo]
            [olivier.math :as om]
            [clojure.math.numeric-tower :as math]))


; Returns the same type it received
(defmulti rev class)

(defmethod rev clojure.lang.PersistentVector [in]
  (vec (rseq in)))

(defmethod rev clojure.lang.PersistentList [in]
  (rseq in))

(defmethod rev :default [in]
  (apply conj () in))

(defn rotations [xs] (take (count xs) (partition (count xs) 1 (cycle xs))))

(defn positions
  "Returns the indexes of all values that match the predicate"
  [pred coll]
  (keep-indexed (fn [idx x] (when (pred x) idx)) coll))

;==============================================================================
;                           Transformations
;
;      from https://groups.google.com/forum/#!topic/clojure/SjmevTjZPcQ
;==============================================================================
(defn conjugate [f g & [g-inv]]
  (comp (or g-inv g) f g))

(defn composite [f f-inv n x]
  (nth (iterate (if (pos? n) f f-inv) x) (math/abs n)))

(defn rotate-left [xs]
  (when (seq xs) (concat (rest xs) [(first xs)])))

(def rotate-right (conjugate rotate-left reverse))

(defn rotate [n xs]
  (composite rotate-left rotate-right n xs))

(defn schwartz [key-fn] #(map (fn [x] [(key-fn x) x]) %))

(def unschwartz #(map second %))

(defn schwartz-sort [key-fn] (conjugate sort (schwartz key-fn) unschwartz))


;;;============================================================================
;;;                           Partitions
;;;============================================================================

(defn partition-by-proportion [x proportion]
  (if (= x 1)
    [1]
    (let [prop (min (max 0 proportion) 1)
          start (min (max 1 (int (* (- 1 prop) x))) (dec x))
          end (- x start)]
      [start end])))



(defn all-subsets
  ([v] (all-subsets v nil))
  ([v len]
   (if (nil? len) (combo/subsets v)
                  (filter #(= (count %) len) (combo/subsets v)))))

(defn all-partitions [n]
  "Returns a vector of all possible integer partitions of n, where n is an integer >= 1"
  (let [parts (combo/partitions (repeat (int n) 1))]
    (vec (filter
           #(> (count %) 1)
           (mapv (fn [v] (mapv #(reduce + %) v)) parts)))))

;;;============================================================================
;;;                           List Transformations
;;;
;;;            (return lists of the same depth as input)
;;;============================================================================

(defn perseveration [v pct advance]
  "Returns a stuttering list"
  (loop [v v coll v]
    (if (empty? v)
      coll
      (let [tail (drop advance v)]
        (recur tail (concat coll (if (< (rand) pct) tail [])))))))

(defn rough-sort
  "Sorts a list into n bins.
   Will produce lists that are mostly ordered, but have some swaps.
   Fewer bins = more swaps."
  ([v bins] (rough-sort v bins <))
  ([v bins sort-fn]
   (let [s (sort-by sort-fn v)
         bins (min (count v) bins)
         full (int (/ (count v) bins))
         parts (mod (count v) bins)
         idxs (sort (flatten (concat (repeat full (range 0 bins)) (take parts (shuffle (range 0 bins))))))
         ; There will be repeated indices in here...
         ]
     (mapv second (sort-by first (shuffle (mapv vector idxs s)))))))


(defn invert
  "Inverts a series of values optionally around a specfied midpoint"
  ([v] (let [min-val (apply min v)
             max-val (apply max v)]
         (map #(- max-val %) (map #(- % min-val) v))))
  ([axis v]
   (map #(- axis (- % axis)) v)))


(defn palindrome
  "Returns a palindromic list/vector depending on type"
  ([coll] (palindrome false coll))
  ([repeat-middle? coll]
   (let [head (if repeat-middle? coll (butlast coll))
         tail (rev coll)]
     (if (vector? coll)
       (into (vec head) tail)
       (into head (rev tail))))))

(defn midpoint-rotate
  "Starts in the middle of a sequence and counts out.  Left?, if specified, determines whether the algorithm starts from
   the left or the right.  Left is the default.

   0 1 2 3 4 5 6 ---> 3 2 4 1 5 0 6 (-> left)
   0 1 2 3 4 5 6 ---> 3 4 2 5 1 6 0 (-> right)
   0 1 2 3 4 5   ---> 3 2 4 1 5 0   (-> left)
   0 1 2 3 4 5   ---> 2 3 1 4 0 5   (-> right)"
  ([v] (midpoint-rotate v true))
  ([v left?]
   (let [len  (count v)
         half (int (/ len 2))
         res  (if (odd? len) [(nth v half)] [])   ;; If odd, treat center value as special case.
         seg1 (reverse (take half v))
         seg2 (drop (- len half) v)
         segs (if left? [seg1 seg2] [seg2 seg1])]
     (into res (apply interleave segs)))))

(defn midpoint-rotations
  "Applies midpoint-rotate to the set until the original would be reached on the next step"
  ([v] (midpoint-rotations v true))
  ([v left?] (take (dec (count v)) (iterate #(midpoint-rotate % left?) v))))

