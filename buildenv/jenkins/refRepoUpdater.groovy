/*******************************************************************************
 * Copyright IBM Corp. and others 2020
 *
 * This program and the accompanying materials are made available under
 * the terms of the Eclipse Public License 2.0 which accompanies this
 * distribution and is available at https://www.eclipse.org/legal/epl-2.0/
 * or the Apache License, Version 2.0 which accompanies this distribution and
 * is available at https://www.apache.org/licenses/LICENSE-2.0.
 *
 * This Source Code may also be made available under the following
 * Secondary Licenses when the conditions for such availability set
 * forth in the Eclipse Public License, v. 2.0 are satisfied: GNU
 * General Public License, version 2 with the GNU Classpath
 * Exception [1] and GNU General Public License, version 2 with the
 * OpenJDK Assembly Exception [2].
 *
 * [1] https://www.gnu.org/software/classpath/license.html
 * [2] http://openjdk.java.net/legal/assembly-exception.html
 *
 * SPDX-License-Identifier: EPL-2.0 OR Apache-2.0 OR GPL-2.0 WITH Classpath-exception-2.0 OR LicenseRef-GPL-2.0 WITH Assembly-exception
 *******************************************************************************/

LABEL = (params.LABEL) ? params.LABEL : '!zOS&&!master&&!proxy'

def jobs = [:]

timeout(time: 6, unit: 'HOURS') {
    timestamps {
        buildNodes = getOnlineNodes()
        for (aNode in buildNodes) {
            def name = aNode
            jobs["${name}"] = {
                node("${name}") {
                    refresh()
                }
            }
        }
        parallel jobs
    }
}

@NonCPS
def getOnlineNodes() {
    def onlineNodes = []
    def buildNodes = jenkins.model.Jenkins.instance.getLabel(LABEL).getNodes()
    for (aNode in buildNodes) {
        if (aNode.toComputer().isOnline()) {
            onlineNodes.add(aNode.getDisplayName())
        }
    }
    return onlineNodes
}

def refresh() {
    dir("${HOME}/gitcache") {
        sh '''
            git init --bare
            git config remote.omr.url https://github.com/eclipse/omr.git
            git config remote.omr.fetch +refs/heads/*:refs/remotes/omr/*
            git fetch omr
        '''
    }
}
