/**
// This is an extension of and edits to TLongLongHashMap from the trove4j library
// downloaded from https://bitbucket.org/trove4j/trove/overview
// 
//   extended to add arrays to hold more values for the same key.
// 
// The original copyright is:
// Copyright (c) 2001, Eric D. Friedman All Rights Reserved.
// Copyright (c) 2009, Rob Eden All Rights Reserved.
// Copyright (c) 2009, Jeff Randall All Rights Reserved.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
// @author Eric D. Friedman
// @author Rob Eden
// @author Jeff Randall
// @version $Id: _K__V_HashMap.template,v 1.1.2.16 2010/03/02 04:09:50 robeden Exp $
*/
package algorithms.util;

import gnu.trove.map.TLongLongMap;
import gnu.trove.procedure.*;
import gnu.trove.set.*;
import gnu.trove.iterator.*;
import gnu.trove.impl.hash.*;
import gnu.trove.impl.HashFunctions;
import gnu.trove.*;
import gnu.trove.set.hash.TLongHashSet;

import java.io.*;
import java.util.*;

/**
 * An open addressed Map implementation for long keys and 
 *    long, long, long, int, int, int values.
 * 
 * Note that the structure was built specifically to hold fields for a tree node,
 * so there are special rules used to help define that a value has not been
 * set.
 * 
 * The key is the unique key for an entry.
 *    value0 is LONG.MIN_VALUE when not used.
 *    value1 is LONG.MIN_VALUE when not used.
 *    value2 is LONG.MIN_VALUE when not used.
 *    value3, as soon as the entry exists, always has a used value.
 *    value4, as soon as the entry exists, always has a used value.
 *    value5, as soon as the entry exists, always has a used value.
 *
 * The fields were built to hold:
 *    value0 is parent key
 *    value1 is left node
 *    value2 is right node
 *    value3 is node value
 *    value4 is color
 *    value5 is size of subtree at key
 * 
 */
public class TLong_LongLongLongIntIntIntHashMap extends TLongLongHash implements 
    Externalizable {
    
    static final long serialVersionUID = 1L;

    /** the values of the map */
    protected transient long[] _values0;
    protected transient long[] _values1;
    protected transient long[] _values2;
    protected transient int[] _values3;
    protected transient int[] _values4;
    protected transient int[] _values5;
    
    protected final long noLinkValue = Long.MIN_VALUE;

    protected final int no_entry_value_int;

    /**
     * Creates a new <code>TLongLongHashMap</code> instance with the default
     * capacity and load factor.
     */
    public TLong_LongLongLongIntIntIntHashMap() {
        super();
        no_entry_value_int = 0;
        no_entry_key = 0;
        no_entry_value = noLinkValue;
    }


    /**
     * Creates a new <code>TLongLongHashMap</code> instance with a prime
     * capacity equal to or greater than <tt>initialCapacity</tt> and
     * with the default load factor.
     *
     * @param initialCapacity an <code>int</code> value
     */
    public TLong_LongLongLongIntIntIntHashMap( int initialCapacity ) {
        super( initialCapacity );
        no_entry_value_int = 0;
        no_entry_key = 0;
        no_entry_value = noLinkValue;
    }


    /**
     * Creates a new <code>TLongLongHashMap</code> instance with a prime
     * capacity equal to or greater than <tt>initialCapacity</tt> and
     * with the specified load factor.
     *
     * @param initialCapacity an <code>int</code> value
     * @param loadFactor a <code>float</code> value
     */
    public TLong_LongLongLongIntIntIntHashMap( int initialCapacity, float loadFactor ) {
        super( initialCapacity, loadFactor );
        no_entry_value_int = 0;
        no_entry_key = 0;
        no_entry_value = noLinkValue;
    }


    /**
     * Creates a new <code>TLongLongHashMap</code> instance with a prime
     * capacity equal to or greater than <tt>initialCapacity</tt> and
     * with the specified load factor.
     * 
     * Note that the "no link value" is fixed to Long.MIN_VALUE.
     *
     * @param initialCapacity an <code>int</code> value
     * @param loadFactor a <code>float</code> value
     * @param noEntryKey a <code>long</code> value that represents
     *                   <tt>null</tt> for the Key set.
     * @param noEntryValue a <code>long</code> value that represents
     *                   <tt>null</tt> for the Value set.
     */
    public TLong_LongLongLongIntIntIntHashMap( int initialCapacity, float loadFactor,
        long noEntryKey, long noEntryValue, int noEntryValueInt ) {
        super( initialCapacity, loadFactor, noEntryKey, Long.MIN_VALUE);
        no_entry_value_int = noEntryValueInt;
        no_entry_key = 0;
        no_entry_value = noLinkValue;
    }


    /**
     * Creates a new <code>TLongLongHashMap</code> instance containing
     * all of the entries in the map passed in.
     *      *
     * @param keys a <tt>long</tt> array containing the keys for the matching values.
     *  
     * @param values0 is <tt>long</tt> array containing the first entry of values.
     *     <em>For trees using this class</em>, this is the <em>parent</em> field.
     *     Note that the backing array has a LONG.MIN_VALUE when the key has
           no parent.
       @param values1 is <tt>long</tt> array containing the second entry of values.
     *     <em>For trees using this class</em>, this is the <em>left</em> field.
     *     Note that the backing array has a LONG.MIN_VALUE when the key has
           no parent.
       @param values2 is <tt>long</tt> array containing the third entry of values.
     *     <em>For trees using this class</em>, this is the <em>right</em> field.
     *     Note that the backing array has a LONG.MIN_VALUE when the key has
           no parent.
       @param values3 is <tt>int</tt> array containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>standard value</em> field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param values4 is <tt>int</tt> array containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>color</em> value field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param values5 <tt>int</tt> array containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>size</em> value field
           holding the size of the subtree of this node.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
     */
    public TLong_LongLongLongIntIntIntHashMap( long[] keys, 
        long[] values0, long[] values1, long[] values2, 
        int[] values3, int[] values4, int[] values5) {
        
        super( Math.max( keys.length, values0.length ) );
        
        no_entry_value_int = 0;
        no_entry_key = 0;
        no_entry_value = noLinkValue;
        
        if (keys.length != values0.length ||
            keys.length != values1.length ||
            keys.length != values2.length ||
            keys.length != values3.length ||
            keys.length != values4.length ||
            keys.length != values5.length) {
            throw new IllegalArgumentException("all arrays must be same length");
        }

        int size = Math.min( keys.length, values0.length );
        for ( int i = 0; i < size; i++ ) {
            this.put( keys[i], 
                values0[i], values1[i], values2[i],
                values3[i], values4[i], values5[i]);
        }
    }
    
    /**
     * from THash, a fast ceiling 
     * @param v
     * @return 
     */
    protected static long fastCeil( double v ) {
        long possible_result = ( long ) v;
        if ( v - possible_result > 0 ) possible_result++;
        return possible_result;
    }

    /* precondition: v > 0 */
    protected static int saturatedCast(long v) {
        int r = (int) (v & 0x7fffffff); // removing sign bit
        if (r != v) {
            return Integer.MAX_VALUE;
        }
        return r;
    }

    /**
     * Creates a new <code>TLongLongHashMap</code> instance containing
     * all of the entries in the map passed in.
     *
     * @param map a <tt>TLongLongMap</tt> that will be duplicated.
     */
    public TLong_LongLongLongIntIntIntHashMap( 
        TLong_LongLongLongIntIntIntHashMap map ) {
        
        super( map.size() );
        
        TLong_LongLongLongIntIntIntHashMap hashmap 
            = ( TLong_LongLongLongIntIntIntHashMap ) map;
        this._loadFactor = Math.abs( hashmap._loadFactor );
        this.no_entry_key = hashmap.no_entry_key;
        this.no_entry_value = hashmap.no_entry_value;
        this.no_entry_value_int = hashmap.no_entry_value_int;
        
        //noinspection RedundantCast
        if ( this.no_entry_key != ( long ) 0 ) {
            Arrays.fill( _set, this.no_entry_key );
        }
                
        setUp( saturatedCast( 
            fastCeil( DEFAULT_CAPACITY / (double) _loadFactor ) ) );
         
        Arrays.fill( _values0, this.noLinkValue );
        Arrays.fill( _values1, this.noLinkValue );
        Arrays.fill( _values2, this.noLinkValue );
        Arrays.fill( _values3, this.no_entry_value_int );
        Arrays.fill( _values4, this.no_entry_value_int );
        Arrays.fill( _values5, this.no_entry_value_int );
       
        putAll( map );
    }


    /**
     * initializes the hashtable to a prime capacity which is at least
     * <tt>initialCapacity + 1</tt>.
     *
     * @param initialCapacity an <code>int</code> value
     * @return the actual capacity chosen
     */
    protected int setUp( int initialCapacity ) {
        int capacity;

        capacity = super.setUp( initialCapacity );
        _values0 = new long[capacity];
        _values1 = new long[capacity];
        _values2 = new long[capacity];
        _values3 = new int[capacity];
        _values4 = new int[capacity];
        _values5 = new int[capacity];
        return capacity;
    }


    /**
     * rehashes the map to the new capacity.
     *
     * @param newCapacity an <code>int</code> value
     */
     /** {@inheritDoc} */
    protected void rehash( int newCapacity ) {
        
        int oldCapacity = _set.length;
        
        long oldKeys[] = _set;
        long oldVals0[] = _values0;
        long oldVals1[] = _values1;
        long oldVals2[] = _values2;
        int oldVals3[] = _values3;
        int oldVals4[] = _values4;
        int oldVals5[] = _values5;
        byte oldStates[] = _states;

        _set = new long[newCapacity];
        _values0 = new long[newCapacity];
        _values1 = new long[newCapacity];
        _values2 = new long[newCapacity];
        _values3 = new int[newCapacity];
        _values4 = new int[newCapacity];
        _values5 = new int[newCapacity];
        _states = new byte[newCapacity];

        for ( int i = oldCapacity; i-- > 0; ) {
            if( oldStates[i] == FULL ) {
                long o = oldKeys[i];
                int index = insertKey( o );
                _values0[index] = oldVals0[i];
                _values1[index] = oldVals1[i];
                _values2[index] = oldVals2[i];
                _values3[index] = oldVals3[i];
                _values4[index] = oldVals4[i];
                _values5[index] = oldVals5[i];
            }
        }
    }

    /**
     * insert key and values.
     * 
     * @param key a <tt>long</tt> containing the key for the associated values.
     *  
     * @param value0 is <tt>long</tt> containing the first entry of values.
     *     <em>For trees using this class</em>, this is the <em>parent</em> field.
     *     Note that the backing array has a LONG.MIN_VALUE when the key has
           no parent.
       @param value1 is <tt>long</tt> containing the second entry of values.
     *     <em>For trees using this class</em>, this is the <em>left</em> field.
     *     Note that the backing array has a LONG.MIN_VALUE when the key has
           no parent.
       @param value2 is <tt>long</tt> containing the third entry of values.
     *     <em>For trees using this class</em>, this is the <em>right</em> field.
     *     Note that the backing array has a LONG.MIN_VALUE when the key has
           no parent.
       @param value3 is <tt>int</tt> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>standard value</em> field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param value4 is <tt>int</tt> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>color</em> value field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param value5 <tt>int</tt> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>size</em> value field
           holding the size of the subtree of this node.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
     */
    public void put( long key, long value0, long value1, long value2,
        int value3, int value4, int value5) {
        int index = insertKey( key );
        doPut( key, value0, value1, value2, value3, value4, value5, index );
    }
    
    /**
     insert key and values.  Note that the missing value0, value1, and value2
     are inserted with the sentinel for long null values, Long.MIN_VALUE.
     * 
     * @param key a <tt>long</tt> containing the key for the associated values.
     *  
       @param value3 is <tt>int</tt> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>standard value</em> field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param value4 is <tt>int</tt> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>color</em> value field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param value5 <tt>int</tt> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>size</em> value field
           holding the size of the subtree of this node.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields. 
     */
    public void put( long key, int value3, int value4, int value5) {
        int index = insertKey( key );
        doPut( key, noLinkValue, noLinkValue, noLinkValue, value3, value4, 
            value5, index );
    }

    /** insert into map if there is not already an existing key.
     * Note that if there is an existing key, false is returned,
     * else true for successful insert.
     * @param key a <tt>long</tt> containing the key for the associated values.
     *  
     * @param value0 is <tt>long</tt> containing the first entry of values.
     *     <em>For trees using this class</em>, this is the <em>parent</em> field.
     *     Note that the backing array has a LONG.MIN_VALUE when the key has
           no parent.
       @param value1 is <tt>long</tt> containing the second entry of values.
     *     <em>For trees using this class</em>, this is the <em>left</em> field.
     *     Note that the backing array has a LONG.MIN_VALUE when the key has
           no parent.
       @param value2 is <tt>long</tt> containing the third entry of values.
     *     <em>For trees using this class</em>, this is the <em>right</em> field.
     *     Note that the backing array has a LONG.MIN_VALUE when the key has
           no parent.
       @param value3 is <tt>int</tt> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>standard value</em> field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param value4 is <tt>int</tt> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>color</em> value field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param value5 <tt>int</tt> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>size</em> value field
           holding the size of the subtree of this node.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
     */
    public boolean putIfAbsent( long key, long value0, long value1, long value2,
        int value3, int value4, int value5 ) {
        int index = insertKey( key );
        if (index < 0)
            return false;
        doPut( key, value0, value1, value2, value3, value4, value5, index );
        return true;
    }
    
    /** insert into map if there is not already an existing key.
     * Note that if there is an existing key, false is returned,
     * else true for successful insert.
     * 
     * Note that the missing value0, value1, and value2
     are inserted with the sentinel for long null values, Long.MIN_VALUE.
     * 
     * @param key a <tt>long</tt> containing the key for the associated values.
     *  
       @param value3 is <tt>int</tt> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>standard value</em> field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param value4 is <tt>int</tt> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>color</em> value field.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields.
       @param value5 <tt>int</tt> containing the fourth entry of values.
     *     <em>For trees using this class</em>, this is the <em>size</em> value field
           holding the size of the subtree of this node.
     *     Note that the backing array should always have a valid value for this
           field, that is, if containsKey(key) is true, this value will be fetched
           and not checked for a null entry like the long fields. 
     */
    public boolean putIfAbsent( long key, int value3, int value4, int value5 ) {
        int index = insertKey( key );
        if (index < 0)
            return false;
        doPut( key, noLinkValue, noLinkValue, noLinkValue, value3, value4, value5, index );
        return true;
    }

    private void doPut( long key, long value0, long value1, long value2,
        int value3, int value4, int value5, int index ) {
        
        boolean isNewMapping = true;
        if ( index < 0 ) {
            index = -index -1;
            isNewMapping = false;
        }
        _values0[index] = value0;
        _values1[index] = value1;
        _values2[index] = value2;
        _values3[index] = value3;
        _values4[index] = value4;
        _values5[index] = value5;

        if (isNewMapping) {
            // a method in THash to expand and rehash if needed
            postInsertHook( consumeFreeSlot );
        }
    }
    
    /**
     * update value0 if if doesn't exist, else return false because there is
     * not an entry already for the given key.
     * @param key
     * @param value0
     * @return 
     */
    public boolean updateValue0( long key, long value0 ) {
        //for existing entries, index returned by insertKey is negative
        int index = insertKey( key );
        if (index >= 0)
            return false;
        if ( index < 0 ) {
            index = -index -1;
        }
        _values0[index] = value0;
        return true;
    }
    /**
     * update value1 if if doesn't exist, else return false because there is
     * not an entry already for the given key.
     * @param key
     * @param value1
     * @return 
     */
    public boolean updateValue1( long key, long value1 ) {
        //for existing entries, index is positive
        int index = index( key );
        if (index < 0) {
            return false;
        }
        _values1[index] = value1;
        return true;
    }
    /**
     * update value2 if if doesn't exist, else return false because there is
     * not an entry already for the given key.
     * @param key
     * @param value2
     * @return 
     */
    public boolean updateValue2( long key, long value2 ) {
        //for existing entries, index is positive
        int index = index( key );
        if (index < 0) {
            return false;
        }
        _values2[index] = value2;
        return true;
    }
    /**
     * update value3 if if doesn't exist, else return false because there is
     * not an entry already for the given key.
     * @param key
     * @param value3
     * @return 
     */
    public boolean updateValue3( long key, int value3 ) {
        //for existing entries, index is positive
        int index = index( key );
        if (index < 0) {
            return false;
        }
        _values3[index] = value3;
        return true;
    }
    /**
     * update value4 if if doesn't exist, else return false because there is
     * not an entry already for the given key.
     * @param key
     * @param value4
     * @return 
     */
    public boolean updateValue4( long key, int value4 ) {
        //for existing entries, index is positive
        int index = index( key );
        if (index < 0) {
            return false;
        }
        _values4[index] = value4;
        return true;
    }
    /**
     * update value5 if if doesn't exist, else return false because there is
     * not an entry already for the given key.
     * @param key
     * @param value5
     * @return 
     */
    public boolean updateValue5( long key, int value5 ) {
        //for existing entries, index is positive
        int index = index( key );
        if (index < 0) {
            return false;
        }
        _values5[index] = value5;
        return true;
    }

    /** {@inheritDoc} */
    public void putAll( TLong_LongLongLongIntIntIntHashMap map ) {
        ensureCapacity( map.size() );
        TLong_LongLongLongIntIntIntHashIterator iter 
            = map.iterator();
        while ( iter.hasNext() ) {
            iter.advance();
            this.put( iter.key(), 
                iter.value0(), iter.value1(), iter.value2(),
                iter.value3(), iter.value4(), iter.value5());
        }
    }


    /** {@inheritDoc} */
    public long getValue0( long key ) {
        int index = index( key );
        return index < 0 ? noLinkValue : _values0[index];
    }
    public long getValue1( long key ) {
        int index = index( key );
        return index < 0 ? noLinkValue : _values1[index];
    }
    public long getValue2( long key ) {
        int index = index( key );
        return index < 0 ? noLinkValue : _values2[index];
    }
    public int getValue3( long key ) {
        int index = index( key );
        return index < 0 ? no_entry_value_int : _values3[index];
    }
    public int getValue4( long key ) {
        int index = index( key );
        return index < 0 ? no_entry_value_int : _values4[index];
    }
    public int getValue5( long key ) {
        int index = index( key );
        return index < 0 ? no_entry_value_int : _values5[index];
    }


    /** {@inheritDoc} */
    public void clear() {
        super.clear();
        int n = _values0.length;
        Arrays.fill( _set, 0, _set.length, no_entry_key );
        Arrays.fill( _values0, 0, n, noLinkValue );
        Arrays.fill( _values1, 0, n, noLinkValue );
        Arrays.fill( _values2, 0, n, noLinkValue );
        Arrays.fill( _values3, 0, n, no_entry_value_int );
        Arrays.fill( _values4, 0, n, no_entry_value_int );
        Arrays.fill( _values5, 0, n, no_entry_value_int );
        Arrays.fill( _states, 0, _states.length, FREE );
    }


    /** {@inheritDoc} */
    public boolean isEmpty() {
        return 0 == _size;
    }


    /** remove the entry for key and return true if succeeded, else false if
     * entry with key was not in map.
     * @param key
     * @return 
     */
    public boolean remove( long key ) {
        long prev = no_entry_value;
        int index = index( key );
        if ( index >= 0 ) {
            prev = _values3[index];
            removeAt( index );    // clear key,state; adjust size
            return true;
        }
        return false;
    }


    /** {@inheritDoc} */
    protected void removeAt( int index ) {
        _values0[index] = noLinkValue;
        _values1[index] = noLinkValue;
        _values2[index] = noLinkValue;
        _values3[index] = no_entry_value_int;
        _values4[index] = no_entry_value_int;
        _values5[index] = no_entry_value_int;
        super.removeAt( index );  // clear key, state; adjust size
    }


    /** {@inheritDoc} */
    public TLongSet keySet() {
        
        TLongSet keys = new TLongHashSet();
        
        if (size() == 0) {
            return keys;
        }
        
        long[] k = _set;
        byte[] states = _states;

        for (int i = k.length, j = 0; i-- > 0;) {
            if (states[i] == FULL) {
                keys.add(k[i]);
            }
        }
        
        return keys;
    }


    /** {@inheritDoc} */
    public long[] keys() {
        long[] keys = new long[size()];
        if ( keys.length == 0 ) {
            return keys;        // nothing to copy
        }
        long[] k = _set;
        byte[] states = _states;

        for ( int i = k.length, j = 0; i-- > 0; ) {
          if ( states[i] == FULL ) {
            keys[j++] = k[i];
          }
        }
        return keys;
    }


    /** {@inheritDoc} */
    public long[] keys( long[] array ) {
        int size = size();
        if ( size == 0 ) {
            return array;       // nothing to copy
        }
        if ( array.length < size ) {
            array = new long[size];
        }

        long[] keys = _set;
        byte[] states = _states;

        for ( int i = keys.length, j = 0; i-- > 0; ) {
          if ( states[i] == FULL ) {
            array[j++] = keys[i];
          }
        }
        return array;
    }

    /** {@inheritDoc} */
    public boolean containsKey( long key ) {
        return contains( key );
    }

    public boolean containsValue0( long key) {
        int index = index( key );
        return (index >= 0);
    }
    public boolean containsValue1( long key) {
        int index = index( key );
        return (index >= 0);
    }
    public boolean containsValue2( long key) {
        int index = index( key );
        return (index >= 0);
    }
    public boolean containsValue3( long key) {
        int index = index( key );
        return (index >= 0);
    }
    public boolean containsValue4( long key) {
        int index = index( key );
        return (index >= 0);
    }
    public boolean containsValue5( long key) {
        int index = index( key );
        return (index >= 0);
    }

    public TLong_LongLongLongIntIntIntHashIterator iterator() {
        return new TLong_LongLongLongIntIntIntHashIterator(this );
    }

    class TLong_LongLongLongIntIntIntHashKeyIterator 
        extends THashPrimitiveIterator implements TLongIterator {

        /**
         * Creates an iterator over the specified map
         *
         * @param hash the <tt>TPrimitiveHash</tt> we will be iterating over.
         */
        TLong_LongLongLongIntIntIntHashKeyIterator( TPrimitiveHash hash ) {
            super( hash );
        }

        /** {@inheritDoc} */
        public long next() {
            moveToNextIndex();
            return _set[_index];
        }

        /** @{inheritDoc} */
        public void remove() {
            if ( _expectedSize != _hash.size() ) {
                throw new ConcurrentModificationException();
            }

            // Disable auto compaction during the remove. This is a workaround for bug 1642768.
            try {
                _hash.tempDisableAutoCompaction();
                TLong_LongLongLongIntIntIntHashMap.this.removeAt( _index );
            }
            finally {
                _hash.reenableAutoCompaction( false );
            }

            _expectedSize--;
        }
    }

    class TLong_LongLongLongIntIntIntHashIterator 
        extends THashPrimitiveIterator {

        /**
         * Creates an iterator over the specified map
         *
         * @param map the <tt>TLongLongHashMap</tt> we will be iterating over.
         */
        TLong_LongLongLongIntIntIntHashIterator( 
            TLong_LongLongLongIntIntIntHashMap map ) {
            super( map );
        }

        /** {@inheritDoc} */
        public void advance() {
            moveToNextIndex();
        }

        /** {@inheritDoc} */
        public long key() {
            return _set[_index];
        }

        /** {@inheritDoc} */
        public long value0() {
            return _values0[_index];
        }
        public long value1() {
            return _values1[_index];
        }
        public long value2() {
            return _values2[_index];
        }
        public int value3() {
            return _values3[_index];
        }
        public int value4() {
            return _values4[_index];
        }
        public int value5() {
            return _values5[_index];
        }

        /** {@inheritDoc} */
        public long setValue0( long val ) {
            long old = value0();
            _values0[_index] = val;
            return old;
        }
        public long setValue1( long val ) {
            long old = value1();
            _values1[_index] = val;
            return old;
        }
        public long setValue2( long val ) {
            long old = value2();
            _values2[_index] = val;
            return old;
        }
        
        public int setValue3( int val ) {
            int old = value3();
            _values3[_index] = val;
            return old;
        }
        public int setValue4( int val ) {
            int old = value4();
            _values4[_index] = val;
            return old;
        }
        public int setValue5( int val ) {
            int old = value5();
            _values5[_index] = val;
            return old;
        }

        /** @{inheritDoc} */
        public void remove() {
            if ( _expectedSize != _hash.size() ) {
                throw new ConcurrentModificationException();
            }
            // Disable auto compaction during the remove. This is a workaround for bug 1642768.
            try {
                _hash.tempDisableAutoCompaction();
                TLong_LongLongLongIntIntIntHashMap.this.removeAt( _index );
            }
            finally {
                _hash.reenableAutoCompaction( false );
            }
            _expectedSize--;
        }
    }


    /** {@inheritDoc}
     * NOTE: this uses just the key and value3 for equals operation.
     */
    @Override
    public boolean equals( Object other ) {
        if ( ! ( other instanceof TLongLongMap ) ) {
            return false;
        }
        TLongLongMap that = ( TLongLongMap ) other;
        if ( that.size() != this.size() ) {
            return false;
        }
        int[] values3 = _values3;
        byte[] states = _states;
        long this_no_entry_value = getNoEntryValue();
        long that_no_entry_value = that.getNoEntryValue();
        for ( int i = values3.length; i-- > 0; ) {
            if ( states[i] == FULL ) {
                long key = _set[i];

                if ( !that.containsKey( key ) ) return false;

                long that_value = that.get( key );
                long this_value = values3[i];
                if ((this_value != that_value)
                    && ( (this_value != this_no_entry_value)
                    || (that_value != that_no_entry_value))
                    ) {

                    return false;
                }
            }
        }
        return true;
    }


    /** {@inheritDoc} 
     * NOTE: this uses the key and values3 for hash code and no other fields.
     */
    @Override
    public int hashCode() {
        int hashcode = 0;
        byte[] states = _states;
        for ( int i = _values0.length; i-- > 0; ) {
            if ( states[i] == FULL ) {
                hashcode += HashFunctions.hash( _set[i] ) ^
                            HashFunctions.hash( _values3[i] );
            }
        }
        return hashcode;
    }


    /** {@inheritDoc} */
    @Override
    public String toString() {

        final StringBuilder buf = new StringBuilder( "{" );
        
        int[] values3 = _values3;
        byte[] states = _states;
        for ( int i = values3.length; i-- > 0; ) {
            if ( states[i] == FULL ) {
                if (buf.length() > 0) {
                    buf.append( ", " );
                }
                long key = _set[i];
                buf.append(key);
                buf.append("=");
                buf.append(_values0[i]);
                buf.append(_values1[i]);
                buf.append(_values2[i]);
                buf.append(_values3[i]);
                buf.append(_values4[i]);
                buf.append(_values5[i]);
            }
        }
        buf.append( "}" );
        return buf.toString();
    }


    /** {@inheritDoc} */
    public void writeExternal(ObjectOutput out) throws IOException {
        // VERSION
    	out.writeByte( 0 );

        // SUPER
    	super.writeExternal( out );

    	// NUMBER OF ENTRIES
    	out.writeInt( _size );

    	// ENTRIES
    	for ( int i = _states.length; i-- > 0; ) {
            if ( _states[i] == FULL ) {
                out.writeLong( _set[i] );
                out.writeLong( _values0[i] );
                out.writeLong( _values1[i] );
                out.writeLong( _values2[i] );
                out.writeInt( _values3[i] );
                out.writeInt( _values4[i] );
                out.writeInt( _values5[i] );
            }
        }
    }


    /** {@inheritDoc} */
    public void readExternal(ObjectInput in) throws IOException, ClassNotFoundException {
        // VERSION
    	in.readByte();

        // SUPER
    	super.readExternal( in );

    	// NUMBER OF ENTRIES
    	int size = in.readInt();
    	setUp( size );

    	// ENTRIES
        while (size-- > 0) {
            long key = in.readLong();
            long val0 = in.readLong();
            long val1 = in.readLong();
            long val2 = in.readLong();
            int val3 = in.readInt();
            int val4 = in.readInt();
            int val5 = in.readInt();
            put(key, val0, val1, val2, val3, val4, val5);
        }
    }
} // TLongLongHashMap
