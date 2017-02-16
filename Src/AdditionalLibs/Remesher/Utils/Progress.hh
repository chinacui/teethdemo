/*===========================================================================*\
 *                                                                           *
 *                              OpenFlipper                                  *
 *           Copyright (c) 2001-2015, RWTH-Aachen University                 *
 *           Department of Computer Graphics and Multimedia                  *
 *                          All rights reserved.                             *
 *                            www.openflipper.org                            *
 *                                                                           *
 *---------------------------------------------------------------------------*
 * This file is part of OpenFlipper.                                         *
 *---------------------------------------------------------------------------*
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  *
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 *                                                                           *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 *                                                                           *
 * 3. Neither the name of the copyright holder nor the names of its          *
 *    contributors may be used to endorse or promote products derived from   *
 *    this software without specific prior written permission.               *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED *
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A           *
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER *
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 *                                                                           *
\*===========================================================================*/

/*
 * Progress.hh
 *
 *  Created on: Mar 14, 2012
 *      Author: ebke
 */

#ifndef PROGRESS_HH_
#define PROGRESS_HH_

#include <set>
#include <algorithm>

#ifdef USE_QT
#include <QProgressDialog>
#endif

namespace ACG {

/**
 * Keep track of the progress of tasks and subtasks.
 *
 * Not thread safe. Could easily be made thread safe, though.
 * See "Thread safety" comments in the source code.
 *
 * Usage Example:
 *
 *     Progress p1(5.0);
 *     p1.increment(1.0);
 *     assert(p1.getNormalizedProgress() == 0.2);
 *     p1.increment(4.0);
 *     assert(p1.getNormalizedProgress() == 1.0);
 *
 *     Progress p2(10.0);
 *     {
 *         Progress p2_1(10.0);
 *         p2.addSubProgress(&p2_1, 5.0);
 *         assert(p2.getNormalizedProgress() == 0.0);
 *         p2_1.increment(5.0);
 *         assert(p2_1.getNormalizedProgress() == 0.5);
 *         assert(p2.getNormalizedProgress() == 0.25);
 *     }
 *
 *     // When the sub progress gets destroyed, the parent
 *     // assumes it is completed.
 *     assert(p2.getNormalizedProgress() == 0.5);
 */

class Progress {
    public:
        Progress(double maxProgress = 1.0) : parent(0), maxProgress(maxProgress), currentProgress(0.0), dialog_padding(0)
        {};
        virtual ~Progress() {
            if (parent) parent->childGoesBye(this);
        };

        void setMaxProgress(double value) {
            maxProgress = value;
            if (currentProgress != 0.0) progressChanged();
        }

#ifdef USE_QT
        void connectToDialog(QProgressDialog *dlg) {
            dialog = dlg;
            dialog->setMaximum(100000);
            updateDialog();
        }
#endif

        /**
         * @return A number in the interval [0, 1], 0 being no progress at all and 1 being a finished job.
         */
        double getNormalizedProgress() const {
            return getProgress() / maxProgress;
        }

        double getProgress() const {
            double result = currentProgress;

            /*
             * Thread safety: Lock mutex on children!
             */
            for (std::set<ChildRecord>::const_iterator it = children.begin(); it != children.end(); ++it) {
                result += it->child->getNormalizedProgress() * it->contribution;
            }

            return std::min(maxProgress, result);
        }

        /**
         * Does NOT take ownership of the supplied Progress object.
         */
        void addSubProgress(Progress *progress, double contribution) {
            progress->parent = this;

            /*
             * Thread safety: Lock mutex on children!
             */
            children.insert(ChildRecord(progress, contribution));
        }

        void increment(double amount) {
            currentProgress += amount;
            progressChanged();
        }

    protected:
        void onChildProgress(Progress *) {
            progressChanged();
        }

        inline void progressChanged() {
            if (parent) parent->onChildProgress(this);
            updateDialog();
        }

        inline void updateDialog() {
#ifdef USE_QT
            if (dialog)
                dialog->setValue(100000 * getNormalizedProgress());
#endif
        }

        inline void childGoesBye(Progress *child) {
            /*
             * Thread safety: Lock mutex on children!
             */
            std::set<ChildRecord>::iterator it = children.find(ChildRecord(child, 0));
            currentProgress += it->contribution;
            children.erase(it);
            progressChanged();
        }

    protected:
        Progress *parent;

        class ChildRecord {
            public:
                ChildRecord(Progress *child, double contribution) : child(child), contribution(contribution) {}

                bool operator==(const ChildRecord &rhs) const { return child == rhs.child; }
                bool operator<(const ChildRecord &rhs) const { return child < rhs.child; }

                Progress *child;
                double contribution;
        };
        std::set<ChildRecord> children;

        double maxProgress, currentProgress;

        union {
#ifdef USE_QT
            QProgressDialog *dialog;
#endif
            /*
             * Make sure the object has the same memory layout, no matter
             * whether USE_QT is defined or not.
             */
            void *dialog_padding;
        };

    private:
        /*
         * Progress objects are not copyable because
         * of the child-parent double linkage.
         */
        Progress(const Progress &rhs);
        Progress &operator=(const Progress &rhs);
};

} /* namespace ACG */
#endif /* PROGRESS_HH_ */
