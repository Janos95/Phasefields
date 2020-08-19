// See permute
template<class BidirIter, class Function>
bool
permute_(BidirIter first1, BidirIter last1,
         typename std::iterator_traits<BidirIter>::difference_type d1,
         Function& f) {
    using std::swap;
    switch(d1) {
        case 0:
        case 1:
            return f();
        case 2:
            if(f())
                return true;
            swap(*first1, *std::next(first1));
            return f();
        case 3: {
            if(f())
                return true;
            BidirIter f2 = std::next(first1);
            BidirIter f3 = std::next(f2);
            swap(*f2, *f3);
            if(f())
                return true;
            swap(*first1, *f3);
            swap(*f2, *f3);
            if(f())
                return true;
            swap(*f2, *f3);
            if(f())
                return true;
            swap(*first1, *f2);
            swap(*f2, *f3);
            if(f())
                return true;
            swap(*f2, *f3);
            return f();
        }
    }
    BidirIter fp1 = std::next(first1);
    for(BidirIter p = fp1; p != last1; ++p){
        if(permute_(fp1, last1, d1 - 1, f))
            return true;
        std::reverse(fp1, last1);
        swap(*first1, *p);
    }
    return permute_(fp1, last1, d1 - 1, f);
}

// Calls f() for each permutation of [first1, last1)
// Divided into permute and permute_ in a (perhaps futile) attempt to
//    squeeze a little more performance out of it.
template<class BidirIter, class Function>
bool
permute(BidirIter first1, BidirIter last1,
        typename std::iterator_traits<BidirIter>::difference_type d1,
        Function& f) {
    using std::swap;
    switch(d1) {
        case 0:
        case 1:
            return f();
        case 2: {
            if(f())
                return true;
            BidirIter i = std::next(first1);
            swap(*first1, *i);
            if(f())
                return true;
            swap(*first1, *i);
        }
            break;
        case 3: {
            if(f())
                return true;
            BidirIter f2 = std::next(first1);
            BidirIter f3 = std::next(f2);
            swap(*f2, *f3);
            if(f())
                return true;
            swap(*first1, *f3);
            swap(*f2, *f3);
            if(f())
                return true;
            swap(*f2, *f3);
            if(f())
                return true;
            swap(*first1, *f2);
            swap(*f2, *f3);
            if(f())
                return true;
            swap(*f2, *f3);
            if(f())
                return true;
            swap(*first1, *f3);
        }
            break;
        default:
            BidirIter fp1 = std::next(first1);
            for(BidirIter p = fp1; p != last1; ++p){
                if(permute_(fp1, last1, d1 - 1, f))
                    return true;
                std::reverse(fp1, last1);
                swap(*first1, *p);
            }
            if(permute_(fp1, last1, d1 - 1, f))
                return true;
            std::reverse(first1, last1);
            break;
    }
    return false;
}
