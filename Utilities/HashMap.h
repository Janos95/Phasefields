

#include <cstdlib>
#include <cstdint>
#include <cstring>

#include <type_traits>

namespace Corrade::Containers {

/* std::pair is not trivially copyable, even if its dependant types are */
template<class K, class V>
struct Element {
    template<class... Args>
    Element(K&& k, Args&& ... args) : key((K&&) k), value((Args&&) args...) {}

    K key;
    V value;
};

template<class T1, class T2>
struct Pair {
    T1 first;
    T2 second;
};

namespace Implementation {

template<class T>
constexpr decltype(auto) move(T&& t) noexcept {
    return static_cast<typename std::remove_reference<T>::type&&>(t);
}

template<class T>
struct IsDestructiveMovable : std::false_type {
};

template<class T> requires std::is_trivial_v<T>
struct IsDestructiveMovable<T> : std::true_type {
};

std::size_t nextPow2(std::uint64_t x) {
    return x == 1ul ? 1ul : 1ul << (64ul - __builtin_clzl(x - 1));
}

template<class K, class T>
struct Bucket {

    //static_assert(std::is_trivial_v<Bucket>);

    using Type = Element<K, T>;

    void swapWithValueInBucket(std::int16_t& dist, std::int16_t& hash, Type& value) {
        using std::swap;
        swap(element(), value);
        swap(distFromIdealBucket, dist);
        //if constexpr(StoresTruncatedHash)
        //    swap(this->hash, hash);
        /* @todo what about end marker */
    }

    template<class... Args>
    void setEmptyBucket(std::int16_t dist, std::int16_t hash, K&& key, Args&& ... args) {
        distFromIdealBucket = dist;
        //if constexpr(StoresTruncatedHash) this->hash = hash;
        ::new(static_cast<void*>(data)) Type{(K&&) key, (Args&&) args...};
    }

    void setEmptyBucket(std::int16_t dist, Bucket& other) noexcept(std::is_nothrow_move_constructible_v<Type>) {
        //if constexpr(StoresTruncatedHash) this->hash = other.hash;
        distFromIdealBucket = dist;
        ::new(static_cast<void*>(data)) Type{(Type&&) other};
    }

    bool empty() { return distFromIdealBucket == EMPTY; }

    void clear() {
        destroy();
        distFromIdealBucket = EMPTY;
    }

    void destroy() { reinterpret_cast<Type*>(data)->~Type(); }

    Type& element() { return reinterpret_cast<Type&>(*data); }

    Type const& element() const { return reinterpret_cast<Type const&>(*data); }

    static constexpr std::int16_t EMPTY = -1;

    alignas(Type) char data[sizeof(Type)];
    std::int16_t distFromIdealBucket;
    bool lastBucket;
};

}

template<class Key, class T, class Hash, int Width = 4>
struct HashMap {

    using BucketType = Implementation::Bucket<Key, T>;
    using Type = typename BucketType::Type;

    template<bool IsConst>
    struct Iterator {

        using ElementType = std::conditional_t<IsConst, const typename BucketType::Type, typename BucketType::Type>;
        BucketType* bucket;

        Iterator& operator++() {
            while(true) {
                if(bucket->lastBucket) {
                    ++bucket;
                    return *this;
                }

                ++bucket;
                if(!bucket->empty()) return *this;
            }
        }

        Iterator operator++(int)& {
            Iterator tmp(*this);
            ++*this;
            return tmp;
        }

        ElementType& operator*() const {
            return bucket->element;
        }

        template<bool IsConst1, bool IsConst2>
        friend bool operator!=(Iterator<IsConst1> const& lhs, Iterator<IsConst2> const& rhs) {
            return lhs.bucket != rhs.bucket;
        }

        template<bool IsConst1, bool IsConst2>
        friend bool operator==(Iterator<IsConst1> const& lhs, Iterator<IsConst2> const& rhs) {
            return lhs.bucket == rhs.bucket;
        }

    };

    using const_iterator = Iterator<true>;
    using iterator = Iterator<false>;

    HashMap() = default;

    explicit HashMap(std::size_t n, Hash const& hash = {}) :
            numBuckets(Implementation::nextPow2(n)), mask(numBuckets - 1), buckets(std::malloc(numBuckets*sizeof(BucketType))) {
        for(int i = 0; i < numBuckets; ++i) {
            buckets[i].distFromIdealBucket = BucketType::Empty;
            buckets[i].lastBucket = false;
        }
    }

    iterator begin() noexcept {
        std::size_t i = 0;
        while(i++ < numBuckets && buckets[i].empty());
        return {buckets + i};
    }

    const_iterator begin() const noexcept { return {const_cast<HashMap*>(this)->begin()}; }

    iterator end() { return {buckets + numBuckets}; }

    const_iterator end() const { return {buckets + numBuckets}; }

    [[nodiscard]] bool empty() const noexcept { return numElements == 0; }

    [[nodiscard]] std::size_t size() const noexcept { return numElements; }

    void eraseFromBucket(iterator pos) {
        pos.bucket->clear();
        numElements--;

        std::size_t previousIdx = pos.bucket - buckets;
        std::size_t idx = nextBucket(previousIdx);

        while(buckets[idx].distFromIdealBucket > 0) {
            const std::int16_t newDistance = buckets[idx].distFromIdealBucket - 1;
            buckets[previousIdx].setEmptyBucket(newDistance, buckets[idx].element());
            buckets[idx].clear();

            previousIdx = idx;
            idx = nextBucket(idx);
        }
    }

    std::size_t nextBucket(std::size_t idx) { return (idx + 1) & mask; }

    template<class... Args>
    Pair<iterator, bool> tryEmplace(Key&& key, Args&& ... args) {
        const std::size_t hash = Hash{}(key);

        std::size_t idx = hash & mask;
        std::int32_t distFromIdealBucket = 0;

        while(distFromIdealBucket <= buckets[idx].distFromIdealBucket) {
            if constexpr(UseStoredHashOnLookUp && StoresTruncatedHash) {
                if(buckets[idx].hash == TruncatedHashType(hash) && key == buckets[idx].key())
                    return {iterator(buckets + idx), false};
            } else if(key == buckets[idx].element().key) {
                return {iterator(buckets + idx), false};
            }

            idx = nextBucket(idx);
            ++distFromIdealBucket;
        }

        if(rehashOnExtremeLoad()) {
            idx = hash & mask;
            distFromIdealBucket = 0;

            while(distFromIdealBucket <= buckets[idx].distFromIdealBucket) {
                idx = nextBucket(idx);
                ++distFromIdealBucket;
            }
        }

        if(buckets[idx].empty()) {
            buckets[idx].setEmptyBucket(distFromIdealBucket, hash, (Key&&) key, (Args&&) args...);
        } else {
            Type value{key, args...};
            buckets[idx].swapWithValueInBucket(distFromIdealBucket, hash, value);
            idx = nextBucket(idx);
            distFromIdealBucket++;

            while(!buckets[idx].empty()) {
                if(distFromIdealBucket > buckets[idx].dist_from_ideal_bucket()) {
                    if(distFromIdealBucket >= maxDistFromIdealBucket)
                        growOnNextInsert = true;

                    buckets[idx].swapWithValueInBucket(distFromIdealBucket, hash, value);
                }

                idx = nextBucket(idx);
                ++distFromIdealBucket;
            }

            buckets[idx].setEmptyBucket(distFromIdealBucket, hash, Implementation::move(value));
        }

        ++numElements;
        return Pair{iterator{buckets + idx}, true};
    }

    T& operator[](Key&& key) { return tryEmplace((Key&&) key).first->value; }

    T const& operator[](Key&& key) const { return tryEmplace((Key&&) key).first->value; }

    iterator find(Key const& key) {


        return end();
    }

    const_iterator find(Key const& key) const { return {const_cast<HashMap*>(this)->find(key).bucket}; }

    static constexpr bool useStoredHashOnRehash(std::size_t numBuckets) {
        if constexpr (StoresTruncatedHash) return numBuckets <= 1 << TruncatedHashType::Width;
        else return false;
    }

    void rehash(std::size_t count) {
        HashMap hashMap(count, Hash{});

        const bool useStoredHash = useStoredHashOnRehash(numBuckets);
        for(int i = 0; i < numBuckets; ++i) {
            if(buckets[i].empty()) continue;

            std::size_t hash;
            if(useStoredHash) hash = buckets[i].hash;
            else hash = Hash{}(buckets[i].key());

            buckets[i].distFromIdealBucket = 0;
            hashMap.insertValueOnRehash(hash & mask, buckets[i]);
        }

        hashMap.numElements = numElements; /* insert value On Rehash does no book keeping */
        hashMap.swap(*this);
        hashMap.destroyOnRehash = true;
    }

    void insertValueOnRehash(std::size_t idx, BucketType& bucket) {
        while(true) {
            if(bucket.distFromIdealBucket > buckets[idx].distFromIdealBucket) {
                if(buckets[idx].empty()) {
                    buckets[idx].set(Implementation::move(bucket));
                    return;
                } else {
                    buckets[idx].swap(bucket);
                }
            }

            ++bucket.distFromIdealBucket;
            idx = nextBucket(idx);
        }
    }

    bool rehashOnExtremeLoad() {
        if(growOnNextInsert || numElements >= loadThreshold) {
            rehash(2*numBuckets);
            growOnNextInsert = false;
            return true;
        }

        if(tryShringOnNextInsert) {
            tryShringOnNextInsert = false;
            if(minLoadFactor != 0.0f && loadFactor() < minLoadFactor) {
                rehash(static_cast<std::size_t>(numElements/maxLoadFactor) + 1);
                return true;
            }
        }
        return false;
    }

    [[nodiscard]] float loadFactor() const {
        if(numBuckets == 0)
            return 0;
        return float(numElements)/float(numBuckets);
    }

    void swap(HashMap& other) {
        using std::swap;
        swap(buckets, other.buckets);
        swap(mask, other.mask);
        swap(numBuckets, other.numBuckets);
        swap(numElements, other.numElements);
        swap(loadThreshold, other.loadThreshold);
        swap(minLoadFactor, other.minLoadFactor);
        swap(maxLoadFactor, other.maxLoadFactor);
        swap(growOnNextInsert, other.growOnNextInsert);
        swap(tryShringOnNextInsert, other.tryShringOnNextInsert);
    }

    ~HashMap() {
        if(!destroyOnRehash) {
            for(int i = 0; i < numBuckets; ++i) {
                if(buckets[i].distFromIdealBucket >= 0)
                    buckets[i].destroy();
            }
        }
        std::free(buckets);
    }

    std::size_t numBuckets;
    std::size_t mask;
    BucketType* buckets;

    std::size_t numElements = 0;

    float minLoadFactor = 0.f;
    float maxLoadFactor = 0.5f;
    std::int16_t maxDistFromIdealBucket = 4096;

    std::size_t loadThreshold;

    bool growOnNextInsert = false;
    bool tryShringOnNextInsert = false;
    bool destroyOnRehash = false;
};

}