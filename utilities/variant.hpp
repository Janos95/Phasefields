//
// Created by janos on 27.04.20.
//

#pragma once

#include <boost/mp11/list.hpp>
#include <boost/mp11/algorithm.hpp>

using namespace boost::mp11;

struct monostate
{
};

constexpr bool operator<(monostate, monostate) noexcept { return false; }
constexpr bool operator>(monostate, monostate) noexcept { return false; }
constexpr bool operator<=(monostate, monostate) noexcept { return true; }
constexpr bool operator>=(monostate, monostate) noexcept { return true; }
constexpr bool operator==(monostate, monostate) noexcept { return true; }
constexpr bool operator!=(monostate, monostate) noexcept { return false; }

template<class... T>
class variant;

template<class D, class... T> union variant_storage_impl;

template<class... T> using variant_storage = variant_storage_impl<mp_all<std::is_trivially_destructible<T>...>, T...>;

template<class D> union variant_storage_impl<D>
{
};

// not all trivially destructible
template<class T1, class... T> union variant_storage_impl<mp_false, T1, T...>
{
    T1 first_;
    variant_storage<T...> rest_;

    template<class... A> constexpr explicit variant_storage_impl( mp_size_t<0>, A&&... a ): first_( std::forward<A>(a)... )
    {
    }

    template<std::size_t I, class... A> constexpr explicit variant_storage_impl( mp_size_t<I>, A&&... a ): rest_( mp_size_t<I-1>(), std::forward<A>(a)... )
    {
    }

    ~variant_storage_impl()
    {
    }

    template<class... A> void emplace( mp_size_t<0>, A&&... a )
    {
        ::new( &first_ ) T1( std::forward<A>(a)... );
    }

    template<std::size_t I, class... A> void emplace( mp_size_t<I>, A&&... a )
    {
        rest_.emplace( mp_size_t<I-1>(), std::forward<A>(a)... );
    }

};

// all trivially destructible
template<class T1, class... T> union variant_storage_impl<mp_true, T1, T...>
{
    T1 first_;
    variant_storage<T...> rest_;

    template<class... A> constexpr explicit variant_storage_impl( mp_size_t<0>, A&&... a ): first_( std::forward<A>(a)... )
    {
    }

    template<std::size_t I, class... A> constexpr explicit variant_storage_impl( mp_size_t<I>, A&&... a ): rest_( mp_size_t<I-1>(), std::forward<A>(a)... )
    {
    }

    template<class... A> void emplace_impl( mp_false, mp_size_t<0>, A&&... a )
    {
        ::new( &first_ ) T1( std::forward<A>(a)... );
    }

    template<std::size_t I, class... A> constexpr void emplace_impl( mp_false, mp_size_t<I>, A&&... a )
    {
        rest_.emplace( mp_size_t<I-1>(), std::forward<A>(a)... );
    }

    template<std::size_t I, class... A> constexpr void emplace_impl( mp_true, mp_size_t<I>, A&&... a )
    {
        *this = variant_storage_impl( mp_size_t<I>(), std::forward<A>(a)... );
    }

    template<std::size_t I, class... A> constexpr void emplace( mp_size_t<I>, A&&... a )
    {
        this->emplace_impl( mp_all<std::is_trivially_move_assignable<T1>, std::is_trivially_move_assignable<T>...>(), mp_size_t<I>(), std::forward<A>(a)... );
    }
};


template<bool is_trivially_destructible, class... T> struct variant_base_impl;

struct none {};

// trivially destructible
template<class... T> struct variant_base_impl<true, T...>
{
    int ix_;
    variant_storage<none, T...> st1_;

    constexpr variant_base_impl(): ix_( 0 ), st1_( mp_size_t<0>() )
    {
    }

    template<class I, class... A> constexpr explicit variant_base_impl( I, A&&... a ): ix_( I::value + 1 ), st1_( mp_size_t<I::value + 1>(), std::forward<A>(a)... )
    {
    }

    // requires: ix_ == 0
    template<class I, class... A> void _replace( I, A&&... a )
    {
        ::new( &st1_ ) variant_storage<none, T...>( mp_size_t<I::value + 1>(), std::forward<A>(a)... );
        ix_ = I::value + 1;
    }

    constexpr std::size_t index() const noexcept
    {
        return ix_ - 1;
    }

    template<std::size_t J, class U, class... A> constexpr void emplace_impl( mp_true, A&&... a )
    {
        static_assert( std::is_nothrow_constructible<U, A&&...>::value, "Logic error: U must be nothrow constructible from A&&..." );

        st1_.emplace( mp_size_t<J>(), std::forward<A>(a)... );
        ix_ = J;
    }

    template<std::size_t J, class U, class... A> constexpr void emplace_impl( mp_false, A&&... a )
    {
        static_assert( std::is_nothrow_move_constructible<U>::value, "Logic error: U must be nothrow move constructible" );

        U tmp( std::forward<A>(a)... );

        st1_.emplace( mp_size_t<J>(), std::move(tmp) );
        ix_ = J;
    }

    template<std::size_t I, class... A> constexpr void emplace( A&&... a )
    {
        std::size_t const J = I+1;
        using U = mp_at_c<variant<T...>, I>;

        this->emplace_impl<J, U>( std::is_nothrow_constructible<U, A&&...>(), std::forward<A>(a)... );
    }
};


// not trivially destructible
template<class... T> struct variant_base_impl<false, T...>
{
    int ix_;
    variant_storage<none, T...> st1_;

    constexpr variant_base_impl(): ix_( 0 ), st1_( mp_size_t<0>() )
    {
    }

    template<class I, class... A> constexpr explicit variant_base_impl( I, A&&... a ): ix_( I::value + 1 ), st1_( mp_size_t<I::value + 1>(), std::forward<A>(a)... )
    {
    }

    // requires: ix_ == 0
    template<class I, class... A> void _replace( I, A&&... a )
    {
        ::new( &st1_ ) variant_storage<none, T...>( mp_size_t<I::value + 1>(), std::forward<A>(a)... );
        ix_ = I::value + 1;
    }

    struct _destroy_L1
    {
        variant_base_impl * this_;

        template<class I> void operator()( I ) const noexcept
        {
            using U = mp_at<mp_list<none, T...>, I>;
            this_->st1_.get( I() ).~U();
        }
    };

    void _destroy() noexcept
    {
        if( ix_ > 0 )
        {
            mp_with_index<1 + sizeof...(T)>( ix_, _destroy_L1{ this } );
        }
    }

    ~variant_base_impl() noexcept
    {
        _destroy();
    }

    constexpr std::size_t index() const noexcept
    {
        return ix_ - 1;
    }

    template<std::size_t I> constexpr mp_at_c<variant<T...>, I>& _get_impl( mp_size_t<I> ) noexcept
    {
        size_t const J = I+1;

        assert( ix_ == J );

        return st1_.get( mp_size_t<J>() );
    }

    template<std::size_t I> constexpr mp_at_c<variant<T...>, I> const& _get_impl( mp_size_t<I> ) const noexcept
    {
        return st1_.get( mp_size_t<I+1>() );
    }

    template<std::size_t I, class... A> void emplace( A&&... a )
    {
        size_t const J = I+1;

        using U = mp_at_c<variant<T...>, I>;

        static_assert( std::is_nothrow_move_constructible<U>::value, "Logic error: U must be nothrow move constructible" );

        U tmp( std::forward<A>(a)... );

        _destroy();

        st1_.emplace( mp_size_t<J>(), std::move(tmp) );
        ix_ = J;
    }
};

