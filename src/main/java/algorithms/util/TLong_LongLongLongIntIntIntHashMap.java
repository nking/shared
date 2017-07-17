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
import gnu.trove.function.TLongFunction;
import gnu.trove.procedure.*;
import gnu.trove.set.*;
import gnu.trove.iterator.*;
import gnu.trove.iterator.hash.*;
import gnu.trove.impl.hash.*;
import gnu.trove.impl.HashFunctions;
import gnu.trove.*;

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
                
        setUp( saturatedCast( fastCeil( DEFAULT_CAPACITY / (double) _loadFactor ) ) );
         
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
    public void putAll( Map<? extends Long, ? extends Long> map ) {
        ensureCapacity( map.size() );
        // could optimize this for cases when map instanceof THashMap
        for ( Map.Entry<? extends Long, ? extends Long> entry : map.entrySet() ) {
            this.put( entry.getKey().longValue(), entry.getValue().longValue() );
        }
    }
    

    /** {@inheritDoc} */
    public void putAll( TLongLongMap map ) {
        ensureCapacity( map.size() );
        TLongLongIterator iter = map.iterator();
        while ( iter.hasNext() ) {
            iter.advance();
            this.put( iter.key(), iter.value() );
        }
    }


    /** {@inheritDoc} */
    public long get( long key ) {
        int index = index( key );
        return index < 0 ? no_entry_value : _values[index];
    }


    /** {@inheritDoc} */
    public void clear() {
        super.clear();
        Arrays.fill( _set, 0, _set.length, no_entry_key );
        Arrays.fill( _values, 0, _values.length, no_entry_value );
        Arrays.fill( _states, 0, _states.length, FREE );
    }


    /** {@inheritDoc} */
    public boolean isEmpty() {
        return 0 == _size;
    }


    /** {@inheritDoc} */
    public long remove( long key ) {
        long prev = no_entry_value;
        int index = index( key );
        if ( index >= 0 ) {
            prev = _values[index];
            removeAt( index );    // clear key,state; adjust size
        }
        return prev;
    }


    /** {@inheritDoc} */
    protected void removeAt( int index ) {
        _values[index] = no_entry_value;
        super.removeAt( index );  // clear key, state; adjust size
    }


    /** {@inheritDoc} */
    public TLongSet keySet() {
        return new TKeyView();
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
    public TLongCollection valueCollection() {
        return new TValueView();
    }


    /** {@inheritDoc} */
    public long[] values() {
        long[] vals = new long[size()];
        if ( vals.length == 0 ) {
            return vals;        // nothing to copy
        }
        long[] v = _values;
        byte[] states = _states;

        for ( int i = v.length, j = 0; i-- > 0; ) {
          if ( states[i] == FULL ) {
            vals[j++] = v[i];
          }
        }
        return vals;
    }


    /** {@inheritDoc} */
    public long[] values( long[] array ) {
        int size = size();
        if ( size == 0 ) {
            return array;       // nothing to copy
        }
        if ( array.length < size ) {
            array = new long[size];
        }

        long[] v = _values;
        byte[] states = _states;

        for ( int i = v.length, j = 0; i-- > 0; ) {
          if ( states[i] == FULL ) {
            array[j++] = v[i];
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

    /** {@inheritDoc} */
    public TLongLongIterator iterator() {
        return new TLongLongHashIterator( this );
    }


    /** {@inheritDoc} */
    public boolean forEachKey( TLongProcedure procedure ) {
        return forEach( procedure );
    }


    /** {@inheritDoc} */
    public boolean forEachValue( TLongProcedure procedure ) {
        byte[] states = _states;
        long[] values = _values;
        for ( int i = values.length; i-- > 0; ) {
            if ( states[i] == FULL && ! procedure.execute( values[i] ) ) {
                return false;
            }
        }
        return true;
    }


    /** {@inheritDoc} */
    public boolean forEachEntry( TLongLongProcedure procedure ) {
        byte[] states = _states;
        long[] keys = _set;
        long[] values = _values;
        for ( int i = keys.length; i-- > 0; ) {
            if ( states[i] == FULL && ! procedure.execute( keys[i], values[i] ) ) {
                return false;
            }
        }
        return true;
    }


    /** {@inheritDoc} */
    public void transformValues( TLongFunction function ) {
        byte[] states = _states;
        long[] values = _values;
        for ( int i = values.length; i-- > 0; ) {
            if ( states[i] == FULL ) {
                values[i] = function.execute( values[i] );
            }
        }
    }


    /** {@inheritDoc} */
    public boolean retainEntries( TLongLongProcedure procedure ) {
        boolean modified = false;
        byte[] states = _states;
        long[] keys = _set;
        long[] values = _values;


        // Temporarily disable compaction. This is a fix for bug #1738760
        tempDisableAutoCompaction();
        try {
            for ( int i = keys.length; i-- > 0; ) {
                if ( states[i] == FULL && ! procedure.execute( keys[i], values[i] ) ) {
                    removeAt( i );
                    modified = true;
                }
            }
        } finally {
            reenableAutoCompaction( true );
        }

        return modified;
    }


    /** {@inheritDoc} */
    public boolean increment( long key ) {
        return adjustValue( key, ( long ) 1 );
    }


    /** {@inheritDoc} */
    public boolean adjustValue( long key, long amount ) {
        int index = index( key );
        if (index < 0) {
            return false;
        } else {
            _values[index] += amount;
            return true;
        }
    }


    /** {@inheritDoc} */
    public long adjustOrPutValue( long key, long adjust_amount, long put_amount ) {
        int index = insertKey( key );
        final boolean isNewMapping;
        final long newValue;
        if ( index < 0 ) {
            index = -index -1;
            newValue = ( _values[index] += adjust_amount );
            isNewMapping = false;
        } else {
            newValue = ( _values[index] = put_amount );
            isNewMapping = true;
        }

        byte previousState = _states[index];

        if ( isNewMapping ) {
            postInsertHook(consumeFreeSlot);
        }

        return newValue;
    }


    /** a view onto the keys of the map. */
    protected class TKeyView implements TLongSet {

        /** {@inheritDoc} */
        public TLongIterator iterator() {
            return new TLongLongKeyHashIterator( TLong_LongLongLongIntIntIntHashMap.this );
        }


        /** {@inheritDoc} */
        public long getNoEntryValue() {
            return no_entry_key;
        }


        /** {@inheritDoc} */
        public int size() {
            return _size;
        }


        /** {@inheritDoc} */
        public boolean isEmpty() {
            return 0 == _size;
        }


        /** {@inheritDoc} */
        public boolean contains( long entry ) {
            return TLong_LongLongLongIntIntIntHashMap.this.contains( entry );
        }


        /** {@inheritDoc} */
        public long[] toArray() {
            return TLong_LongLongLongIntIntIntHashMap.this.keys();
        }


        /** {@inheritDoc} */
        public long[] toArray( long[] dest ) {
            return TLong_LongLongLongIntIntIntHashMap.this.keys( dest );
        }


        /**
         * Unsupported when operating upon a Key Set view of a TLongLongMap
         * <p/>
         * {@inheritDoc}
         */
        public boolean add( long entry ) {
            throw new UnsupportedOperationException();
        }


        /** {@inheritDoc} */
        public boolean remove( long entry ) {
            return no_entry_value != TLong_LongLongLongIntIntIntHashMap.this.remove( entry );
        }


        /** {@inheritDoc} */
        public boolean containsAll( Collection<?> collection ) {
            for ( Object element : collection ) {
                if ( element instanceof Long ) {
                    long ele = ( ( Long ) element ).longValue();
                    if ( ! TLong_LongLongLongIntIntIntHashMap.this.containsKey( ele ) ) {
                        return false;
                    }
                } else {
                    return false;
                }
            }
            return true;
        }


        /** {@inheritDoc} */
        public boolean containsAll( TLongCollection collection ) {
            TLongIterator iter = collection.iterator();
            while ( iter.hasNext() ) {
                if ( ! TLong_LongLongLongIntIntIntHashMap.this.containsKey( iter.next() ) ) {
                    return false;
                }
            }
            return true;
        }


        /** {@inheritDoc} */
        public boolean containsAll( long[] array ) {
            for ( long element : array ) {
                if ( ! TLong_LongLongLongIntIntIntHashMap.this.contains( element ) ) {
                    return false;
                }
            }
            return true;
        }


        /**
         * Unsupported when operating upon a Key Set view of a TLongLongMap
         * <p/>
         * {@inheritDoc}
         */
        public boolean addAll( Collection<? extends Long> collection ) {
            throw new UnsupportedOperationException();
        }


        /**
         * Unsupported when operating upon a Key Set view of a TLongLongMap
         * <p/>
         * {@inheritDoc}
         */
        public boolean addAll( TLongCollection collection ) {
            throw new UnsupportedOperationException();
        }


        /**
         * Unsupported when operating upon a Key Set view of a TLongLongMap
         * <p/>
         * {@inheritDoc}
         */
        public boolean addAll( long[] array ) {
            throw new UnsupportedOperationException();
        }


        /** {@inheritDoc} */
        @SuppressWarnings({"SuspiciousMethodCalls"})
        public boolean retainAll( Collection<?> collection ) {
            boolean modified = false;
            TLongIterator iter = iterator();
            while ( iter.hasNext() ) {
                if ( ! collection.contains( Long.valueOf ( iter.next() ) ) ) {
                    iter.remove();
                    modified = true;
                }
            }
            return modified;
        }


        /** {@inheritDoc} */
        public boolean retainAll( TLongCollection collection ) {
            if ( this == collection ) {
                return false;
            }
            boolean modified = false;
            TLongIterator iter = iterator();
            while ( iter.hasNext() ) {
                if ( ! collection.contains( iter.next() ) ) {
                    iter.remove();
                    modified = true;
                }
            }
            return modified;
        }


        /** {@inheritDoc} */
        public boolean retainAll( long[] array ) {
            boolean changed = false;
            Arrays.sort( array );
            long[] set = _set;
            byte[] states = _states;

            for ( int i = set.length; i-- > 0; ) {
                if ( states[i] == FULL && ( Arrays.binarySearch( array, set[i] ) < 0) ) {
                    removeAt( i );
                    changed = true;
                }
            }
            return changed;
        }


        /** {@inheritDoc} */
        public boolean removeAll( Collection<?> collection ) {
            boolean changed = false;
            for ( Object element : collection ) {
                if ( element instanceof Long ) {
                    long c = ( ( Long ) element ).longValue();
                    if ( remove( c ) ) {
                        changed = true;
                    }
                }
            }
            return changed;
        }


        /** {@inheritDoc} */
        public boolean removeAll( TLongCollection collection ) {
            if ( this == collection ) {
                clear();
                return true;
            }
            boolean changed = false;
            TLongIterator iter = collection.iterator();
            while ( iter.hasNext() ) {
                long element = iter.next();
                if ( remove( element ) ) {
                    changed = true;
                }
            }
            return changed;
        }


        /** {@inheritDoc} */
        public boolean removeAll( long[] array ) {
            boolean changed = false;
            for ( int i = array.length; i-- > 0; ) {
                if ( remove( array[i] ) ) {
                    changed = true;
                }
            }
            return changed;
        }


        /** {@inheritDoc} */
        public void clear() {
            TLong_LongLongLongIntIntIntHashMap.this.clear();
        }


        /** {@inheritDoc} */
        public boolean forEach( TLongProcedure procedure ) {
            return TLong_LongLongLongIntIntIntHashMap.this.forEachKey( procedure );
        }


        @Override
        public boolean equals( Object other ) {
            if (! (other instanceof TLongSet)) {
                return false;
            }
            final TLongSet that = ( TLongSet ) other;
            if ( that.size() != this.size() ) {
                return false;
            }
            for ( int i = _states.length; i-- > 0; ) {
                if ( _states[i] == FULL ) {
                    if ( ! that.contains( _set[i] ) ) {
                        return false;
                    }
                }
            }
            return true;
        }


        @Override
        public int hashCode() {
            int hashcode = 0;
            for ( int i = _states.length; i-- > 0; ) {
                if ( _states[i] == FULL ) {
                    hashcode += HashFunctions.hash( _set[i] );
                }
            }
            return hashcode;
        }


        @Override
        public String toString() {
            final StringBuilder buf = new StringBuilder( "{" );
            forEachKey( new TLongProcedure() {
                private boolean first = true;


                public boolean execute( long key ) {
                    if ( first ) {
                        first = false;
                    } else {
                        buf.append( ", " );
                    }

                    buf.append( key );
                    return true;
                }
            } );
            buf.append( "}" );
            return buf.toString();
        }
    }


    /** a view onto the values of the map. */
    protected class TValueView implements TLongCollection {

        /** {@inheritDoc} */
        public TLongIterator iterator() {
            return new TLongLongValueHashIterator( TLong_LongLongLongIntIntIntHashMap.this );
        }


        /** {@inheritDoc} */
        public long getNoEntryValue() {
            return no_entry_value;
        }


        /** {@inheritDoc} */
        public int size() {
            return _size;
        }


        /** {@inheritDoc} */
        public boolean isEmpty() {
            return 0 == _size;
        }


        /** {@inheritDoc} */
        public boolean contains( long entry ) {
            return TLong_LongLongLongIntIntIntHashMap.this.containsValue( entry );
        }


        /** {@inheritDoc} */
        public long[] toArray() {
            return TLong_LongLongLongIntIntIntHashMap.this.values();
        }


        /** {@inheritDoc} */
        public long[] toArray( long[] dest ) {
            return TLong_LongLongLongIntIntIntHashMap.this.values( dest );
        }



        public boolean add( long entry ) {
            throw new UnsupportedOperationException();
        }


        /** {@inheritDoc} */
        public boolean remove( long entry ) {
            long[] values = _values;
            byte[] states = _states;

            for ( int i = values.length; i-- > 0; ) {
                if ( ( states[i] != FREE && states[i] != REMOVED ) && entry == values[i] ) {
                    removeAt( i );
                    return true;
                }
            }
            return false;
        }

        /** {@inheritDoc} */
        public boolean containsAll( Collection<?> collection ) {
            for ( Object element : collection ) {
                if ( element instanceof Long ) {
                    long ele = ( ( Long ) element ).longValue();
                    if ( ! TLong_LongLongLongIntIntIntHashMap.this.containsValue( ele ) ) {
                        return false;
                    }
                } else {
                    return false;
                }
            }
            return true;
        }


        /** {@inheritDoc} */
        public boolean containsAll( TLongCollection collection ) {
            TLongIterator iter = collection.iterator();
            while ( iter.hasNext() ) {
                if ( ! TLong_LongLongLongIntIntIntHashMap.this.containsValue( iter.next() ) ) {
                    return false;
                }
            }
            return true;
        }


        /** {@inheritDoc} */
        public boolean containsAll( long[] array ) {
            for ( long element : array ) {
                if ( ! TLong_LongLongLongIntIntIntHashMap.this.containsValue( element ) ) {
                    return false;
                }
            }
            return true;
        }


        /** {@inheritDoc} */
        public boolean addAll( Collection<? extends Long> collection ) {
            throw new UnsupportedOperationException();
        }


        /** {@inheritDoc} */
        public boolean addAll( TLongCollection collection ) {
            throw new UnsupportedOperationException();
        }


        /** {@inheritDoc} */
        public boolean addAll( long[] array ) {
            throw new UnsupportedOperationException();
        }


        /** {@inheritDoc} */
        @SuppressWarnings({"SuspiciousMethodCalls"})
        public boolean retainAll( Collection<?> collection ) {
            boolean modified = false;
            TLongIterator iter = iterator();
            while ( iter.hasNext() ) {
                if ( ! collection.contains( Long.valueOf ( iter.next() ) ) ) {
                    iter.remove();
                    modified = true;
                }
            }
            return modified;
        }


        /** {@inheritDoc} */
        public boolean retainAll( TLongCollection collection ) {
            if ( this == collection ) {
                return false;
            }
            boolean modified = false;
            TLongIterator iter = iterator();
            while ( iter.hasNext() ) {
                if ( ! collection.contains( iter.next() ) ) {
                    iter.remove();
                    modified = true;
                }
            }
            return modified;
        }


        /** {@inheritDoc} */
        public boolean retainAll( long[] array ) {
            boolean changed = false;
            Arrays.sort( array );
            long[] values = _values;
            byte[] states = _states;

            for ( int i = values.length; i-- > 0; ) {
                if ( states[i] == FULL && ( Arrays.binarySearch( array, values[i] ) < 0) ) {
                    removeAt( i );
                    changed = true;
                }
            }
            return changed;
        }


        /** {@inheritDoc} */
        public boolean removeAll( Collection<?> collection ) {
            boolean changed = false;
            for ( Object element : collection ) {
                if ( element instanceof Long ) {
                    long c = ( ( Long ) element ).longValue();
                    if ( remove( c ) ) {
                        changed = true;
                    }
                }
            }
            return changed;
        }


        /** {@inheritDoc} */
        public boolean removeAll( TLongCollection collection ) {
            if ( this == collection ) {
                clear();
                return true;
            }
            boolean changed = false;
            TLongIterator iter = collection.iterator();
            while ( iter.hasNext() ) {
                long element = iter.next();
                if ( remove( element ) ) {
                    changed = true;
                }
            }
            return changed;
        }


        /** {@inheritDoc} */
        public boolean removeAll( long[] array ) {
            boolean changed = false;
            for ( int i = array.length; i-- > 0; ) {
                if ( remove( array[i] ) ) {
                    changed = true;
                }
            }
            return changed;
        }


        /** {@inheritDoc} */
        public void clear() {
            TLong_LongLongLongIntIntIntHashMap.this.clear();
        }


        /** {@inheritDoc} */
        public boolean forEach( TLongProcedure procedure ) {
            return TLong_LongLongLongIntIntIntHashMap.this.forEachValue( procedure );
        }


        /** {@inheritDoc} */
        @Override
        public String toString() {
            final StringBuilder buf = new StringBuilder( "{" );
            forEachValue( new TLongProcedure() {
                private boolean first = true;

                public boolean execute( long value ) {
                    if ( first ) {
                        first = false;
                    } else {
                        buf.append( ", " );
                    }

                    buf.append( value );
                    return true;
                }
            } );
            buf.append( "}" );
            return buf.toString();
        }
    }


    class TLongLongKeyHashIterator extends THashPrimitiveIterator implements TLongIterator {

        /**
         * Creates an iterator over the specified map
         *
         * @param hash the <tt>TPrimitiveHash</tt> we will be iterating over.
         */
        TLongLongKeyHashIterator( TPrimitiveHash hash ) {
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


   
    class TLongLongValueHashIterator extends THashPrimitiveIterator implements TLongIterator {

        /**
         * Creates an iterator over the specified map
         *
         * @param hash the <tt>TPrimitiveHash</tt> we will be iterating over.
         */
        TLongLongValueHashIterator( TPrimitiveHash hash ) {
            super( hash );
        }

        /** {@inheritDoc} */
        public long next() {
            moveToNextIndex();
            return _values[_index];
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


    class TLongLongHashIterator extends THashPrimitiveIterator implements TLongLongIterator {

        /**
         * Creates an iterator over the specified map
         *
         * @param map the <tt>TLongLongHashMap</tt> we will be iterating over.
         */
        TLongLongHashIterator( TLong_LongLongLongIntIntIntHashMap map ) {
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
        public long value() {
            return _values[_index];
        }

        /** {@inheritDoc} */
        public long setValue( long val ) {
            long old = value();
            _values[_index] = val;
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


    /** {@inheritDoc} */
    @Override
    public boolean equals( Object other ) {
        if ( ! ( other instanceof TLongLongMap ) ) {
            return false;
        }
        TLongLongMap that = ( TLongLongMap ) other;
        if ( that.size() != this.size() ) {
            return false;
        }
        long[] values = _values;
        byte[] states = _states;
        long this_no_entry_value = getNoEntryValue();
        long that_no_entry_value = that.getNoEntryValue();
        for ( int i = values.length; i-- > 0; ) {
            if ( states[i] == FULL ) {
                long key = _set[i];

                if ( !that.containsKey( key ) ) return false;

                long that_value = that.get( key );
                long this_value = values[i];
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


    /** {@inheritDoc} */
    @Override
    public int hashCode() {
        int hashcode = 0;
        byte[] states = _states;
        for ( int i = _values.length; i-- > 0; ) {
            if ( states[i] == FULL ) {
                hashcode += HashFunctions.hash( _set[i] ) ^
                            HashFunctions.hash( _values[i] );
            }
        }
        return hashcode;
    }


    /** {@inheritDoc} */
    @Override
    public String toString() {
        final StringBuilder buf = new StringBuilder( "{" );
        forEachEntry( new TLongLongProcedure() {
            private boolean first = true;
            public boolean execute( long key, long value ) {
                if ( first ) first = false;
                else buf.append( ", " );

                buf.append(key);
                buf.append("=");
                buf.append(value);
                return true;
            }
        });
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
                out.writeLong( _values[i] );
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
            long val = in.readLong();
            put(key, val);
        }
    }
} // TLongLongHashMap
