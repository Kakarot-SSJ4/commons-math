/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package org.apache.commons.math4.stat.descriptive;

import org.apache.commons.math4.exception.MathIllegalArgumentException;
import org.apache.commons.math4.util.MathArrays;

import org.checkerframework.checker.index.qual.NonNegative;
import org.checkerframework.checker.index.qual.IndexOrHigh;
import org.checkerframework.checker.index.qual.LTLengthOf;
import org.checkerframework.common.value.qual.MinLen;


/**
 * Base interface implemented by all statistics.
 */
public interface UnivariateStatistic extends MathArrays.Function {
    /**
     * Returns the result of evaluating the statistic over the input array.
     *
     * @param values input array
     * @return the value of the statistic applied to the input array
     * @throws MathIllegalArgumentException  if values is null
     */
    @Override
    double evaluate(double[] values) throws MathIllegalArgumentException;

    /**
     * Returns the result of evaluating the statistic over the specified entries
     * in the input array.
     *
     * @param values the input array
     * @param begin the index of the first element to include
     * @param length the number of elements to include
     * @return the value of the statistic applied to the included array entries
     * @throws MathIllegalArgumentException if values is null or the indices are invalid
     */
    @Override
    double evaluate(double[] values, @IndexOrHigh("#1") int begin, @NonNegative @LTLengthOf(value = {"#1"}, offset = {"#2 - 1"}) int length) throws MathIllegalArgumentException;

    /**
     * Returns a copy of the statistic with the same internal state.
     *
     * @return a copy of the statistic
     */
    UnivariateStatistic copy();
}
